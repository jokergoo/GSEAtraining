




tb = read.table("E-GEOD-101794-analytics.tsv", sep = "\t", skip = 1)

df = data.frame(gene = tb[, 2], fdr = p.adjust(tb[, 3], "BH"))

l = df$fdr < 0.05
l[is.na(l)] = FALSE

diff_gene = df$gene[l]
diff_gene = diff_gene[diff_gene != ""]
diff_gene = unique(diff_gene)


bg_gene = df$gene[df$gene != ""]
bg_gene = unique(bg_gene)

ln = readLines("h.all.v7.5.1.symbols.gmt")
ln = strsplit(ln, "\t")

gs = lapply(ln, function(x) x[-(1:2)])
names(gs) = lapply(ln, function(x) x[1])


fisher = function(x, y, bg) {
	n11 = length(intersect(x, y))
	n = length(bg)

	n10 = length(x)
	n01 = length(y)

	n12 = n10 - n11
	n21 = n01 - n11
	n22 = n - 10 - n21

	fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))$p.value
}



gene_set = gs[["HALLMARK_KRAS_SIGNALING_DN"]]



1 - phyper(13, 200, 38392, 968)
1 - phyper(13, 968, 37624, 200)


fisher.test(matrix(c(14, 186, 954, 37438), nrow = 2))

1 - pbinom(13, 200, 968/38592)
1 - pbinom(13, 968, 200/38592)

chisq.test(matrix(c(14, 186, 954, 37438), nrow = 2), correct = FALSE)

p1 = 14/200
p2 = 954/38392
p = 968/38592

z = abs(p1 - p2)/sqrt(p*(1-p))/sqrt(1/200 + 1/38392)
2*pnorm(z, lower.tail = FALSE)
z^2


p1 = 14/968
p2 = 186/37624
p = 200/38592

z = abs(p1 - p2)/sqrt(p*(1-p))/sqrt(1/968 + 1/37624)
2*pnorm(z, lower.tail = FALSE)
z^2





library(microbenchmark)

microbenchmark(
	hyper = 1 - phyper(13, 200, 38392, 968),
	fisher = fisher.test(matrix(c(14, 186, 954, 37438), nrow = 2)),
	binom = 1 - pbinom(13, 968, 200/38592),
	chisq = chisq.test(matrix(c(14, 186, 954, 37438), nrow = 2), correct = FALSE),
	ztest = {
		p1 = 14/200
		p2 = 954/38392
		p = 968/38592

		z = abs(p1 - p2)/sqrt(p*(1-p))/sqrt(1/200 + 1/38392)
		2*pnorm(z, lower.tail = FALSE)
	},
	times = 1000
)



1 - phyper(13, 200, 38392, 968)
1 - phyper(13, 200, 38392 - 10000, 968)
1 - phyper(13, 200, 38392 - 20000, 968)
1 - phyper(13, 200, 38392 + 10000, 968)
1 - phyper(13, 200, 38392 + 20000, 968)


1 - phyper(13, 200, 38392, 968)
1 - phyper(13-3, 200, 38392, 968)
1 - phyper(13+3, 200, 38392, 968)



####### GSEA ######

condition = read.cls("P53.cls", treatment = "MUT", control = "WT")$label
expr = read.gct("P53_collapsed_symbols.gct")

ln = strsplit(readLines("c2.symbols.gmt"), "\t")
gs = lapply(ln, function(x) x[-(1:2)])
names(gs) = sapply(ln, function(x) x[1])

geneset = gs[["p53hypoxiaPathway"]]

s = apply(expr, 1, function(x) {
	x1 = x[condition == "WT"]
	x2 = x[condition == "MUT"]
	(mean(x1) - mean(x2))/(sd(x2) + sd(x1))
})

s = sort(s, decreasing = TRUE)

## original GSEA
l_set = names(s) %in% geneset
f1 = cumsum(l_set)/sum(l_set)

l_other = !names(s) %in% geneset
f2 = cumsum(l_other)/sum(l_other)

n = length(s)

plot(1:n, f1, type = "l", col = "red")
lines(1:n, f2, col = "blue")

plot(f1 - f2, type = "l")
abline(h = 0, lty = 2, col = "grey")
points(which(l_set), rep(0, sum(l_set)), pch = "|", col = "red")

es = max(f1 - f2)

## directly test via KS-test

ks.test(which(l_set), which(l_other))


## by sample permutation
library(matrixStats)
calculate_es = function(condition) {

	m1 = expr[, condition == "WT"]
	m2 = expr[, condition == "MUT"]

	s = (rowMeans(m1) - rowMeans(m2))/(rowSds(m1) + rowSds(m2))

	s = sort(s, decreasing = TRUE)

	l_set = names(s) %in% geneset
	f1 = cumsum(l_set)/sum(l_set)

	l_other = !l_set
	f2 = cumsum(l_other)/sum(l_other)

	max(f1 - f2)
}


es = calculate_es(condition)

es_rand = numeric(1000)
for(i in 1:1000) {
	es_rand[i] = calculate_es(sample(condition, length(condition)))
}
sum(es_rand >= es)/1000


hist(es_rand)
abline(v = es, col = "red")

### improved version

calculate_es_v2 = function(condition, plot = FALSE, power = 1) {

	m1 = expr[, condition == "WT"]
	m2 = expr[, condition == "MUT"]

	s = (rowMeans(m1) - rowMeans(m2))/(rowSds(m1) + rowSds(m2))

	s = sort(s, decreasing = TRUE)

	l_set = names(s) %in% geneset
	s_set = abs(s)^power
	s_set[!l_set] = 0
	f1 = cumsum(s_set)/sum(s_set)

	l_other = !l_set
	f2 = cumsum(l_other)/sum(l_other)

	if(plot) {
		plot(f1 - f2, type = "l")
		abline(h = 0, lty = 2, col = "grey")
		points(which(l_set), rep(0, sum(l_set)), pch = "|", col = "red")
	}

	max(f1 - f2)
}

es = calculate_es_v2(condition)
es_rand = numeric(1000)
for(i in 1:1000) {
	es_rand[i] = calculate_es_v2(sample(condition, length(condition)))
}

sum(es_rand >= es)/1000


hist(es_rand, xlim = c(0, 1))
abline(v = es, col = "red")

### by gene permutation


calculate_es_v2_gene_perm = function(s, perm = FALSE, power = 1) {

	s = sort(s, decreasing = TRUE)
	if(perm) {
		names(s) = sample(names(s), length(s))
	}

	l_set = names(s) %in% geneset
	s_set = abs(s)^power
	s_set[!l_set] = 0
	f1 = cumsum(s_set)/sum(s_set)

	l_other = !l_set
	f2 = cumsum(l_other)/sum(l_other)

	max(f1 - f2)
}


m1 = expr[, condition == "WT"]
m2 = expr[, condition == "MUT"]

s = (rowMeans(m1) - rowMeans(m2))/(rowSds(m1) + rowSds(m2))

es = calculate_es_v2_gene_perm(s)
es_rand = numeric(1000)
for(i in 1:1000) {
	es_rand[i] = calculate_es_v2_gene_perm(s, perm = TRUE)
}

sum(es_rand >= es)/1000


hist(es_rand, xlim = c(0, 1))
abline(v = es, col = "red")


## clusterProfiler
cola:::guess_id_mapping(names(s)) -> map
s2 = s
names(s2) = map[names(s)]
s2 = s2[!is.na(names(s2))]
gene_diff_score = s2
gene_diff_score = sort(gene_diff_score, decre = T)
res = gseGO(geneList=gene_diff_score, OrgDb = org.Hs.eg.db, ont = "BP")
gseaplot2(res, geneSetID=1)



convert_to_entrez_id = function(x) {
	if(is.matrix(x)) {
		map = cola:::guess_id_mapping(rownames(x))
		new_rn = map[rownames(x)]
		l = is.na(new_rn)

		x = x[!l, , drop = FALSE]
		new_rn = new_rn[!l]

		x2 = do.call(rbind, tapply(1:nrow(x), new_rn, function(ind) {
			colMeans(x[ind, , drop = FALSE])
		}))
		return(x2)

	} else {
		map = cola:::guess_id_mapping(names(x))
		x2 = s
		names(x2) = map[names(x)]
		x2 = x2[!is.na(names(x2))]
		x2 = tapply(x2, names(x2), mean)
		return(x2)
	}
}
