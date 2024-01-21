

#' Convert gene sets between list and a data frame
#' 
#' @param list A list of genes.
#' @rdname list_to_dataframe
#' @export
list_to_dataframe = function(list) {
    df = data.frame(gene_set = rep(names(list), times = sapply(list, length)),
               gene = unlist(list))
    rownames(df) = NULL
    df
}


#' @param df A two-column data frame
#' 
#' @details
#' In `gs_dataframe_to_list()`, which column contains genes and which column contains gene sets
#' are automatically checked by the number of genes and gene sets. Basically number of genes
#' should be larger than the number of gene sets.
#' 
#' @rdname list_to_dataframe
#' @export
dataframe_to_list = function(df) {
    n1 = length(unique(df[, 1]))
    n2 = length(unique(df[, 2]))
    if(n1 < n2) {
        split(df[, 2], df[, 1])
    } else {
        split(df[, 1], df[, 2])
    }
}

msigdb_base_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
msigdb_env = new.env(parent = emptyenv())
msigdb_env$all_versions = NULL
msigdb_env$file_list = list()
msigdb_env$all_collections = list()
msigdb_env$gene_sets = list()


#' Retrieve gene sets from MSigDB
#' 
#' @rdname msigdb
#' @export
list_msigdb_versions = function() {
    if(is.null(msigdb_env$all_versions)) {
        all_versions = get_file_list(msigdb_base_url)
        msigdb_env$all_versions = all_versions
    } else {
        all_versions = msigdb_env$all_versions
    }
    all_versions
}

#' @importFrom utils menu
choose_msigdb_version = function() {
    all_versions = list_msigdb_versions()
    ind = menu(all_versions, title = "Choose an MSigDB version:")
    all_versions[ind]
}

#' @rdname msigdb
#' 
#' @param version The MSigDB version. The value can be obtained by `list_msigdb_versions()`. 
#'         If this argument is not specified, a menu will be printed to let users to select.
#' @export
list_msigdb_collections = function(version = NULL) {
    if(is.null(version)) {
        version = choose_msigdb_version()
    }

    if(is.null(msigdb_env$file_list[[version]])) {
        files = get_file_list(paste0(msigdb_base_url, "/", version))
        files = files[grep("\\.gmt$", files)]
        files = files[!grepl("^msigdb", files)]
        msigdb_env$file_list[[version]] = files
    } else {
        files = msigdb_env$file_list[[version]]
    }
    collections = unique(gsub(paste0(".v", version, ".*$"), "", files))
    msigdb_env$all_collections[[version]] = collections
    collections
}

choose_msigdb_collection = function(version) {
    
    collections = list_msigdb_collections(version)
    ind = menu(collections, title = paste0("Choose a gene set collection for version ", version, ":"))
    collections[ind]
}

get_file_list = function(url) {
    con = url(url)
    on.exit(close(con))
    ln = readLines(con)
    ind = grep("^<tr><td", ln)
    rows = ln[ind]
    rows = gsub("</td><td[^>]*>", ";", rows)
    rows = gsub("<.*?>", "", rows)
    files = sapply(strsplit(rows, ";"), function(x) x[2])[-1]
    gsub("/", "", files)
}


#' @rdname msigdb
#' 
#' @param collection The gene set collection. The values can be obtained by `list_msigdb_collections()`.
#' @param gene_id_type One of "entrez" and "symbols".
#' @param as_table Whether to return a list or a table.
#' @export
get_msigdb = function(version = choose_msigdb_version(), 
    collection = choose_msigdb_collection(version),
    gene_id_type = c("entrez", "symbols"), as_table = FALSE) {
    
    version = force(version)
    gene_id_type = match.arg(gene_id_type)

    if(is.null(msigdb_env$all_versions)) {
        list_msigdb_versions()
    }
    
    if(!version %in% msigdb_env$all_versions) {
        i = grep(version, msigdb_env$all_versions, ignore.case = TRUE)
        if(length(i) == 1) {
            version = msigdb_env$all_versions[i]
        } else {
            message(paste0("Cannot find version '", version, "', please select a valid value."))
            version = choose_msigdb_version()
        }
    }

    collection = force(collection)
    if(is.null(msigdb_env$all_collections[[version]])) {
        list_msigdb_collections(version)
    }
    if(!collection %in% msigdb_env$all_collections[[version]]) {
        i = grep(collection, msigdb_env$all_collections[[version]], ignore.case = TRUE)
        if(length(i) == 1) {
            collection = msigdb_env$all_collections[[version]][i]
        } else {
            message(paste0("Cannot find collection '", collection, "', please select a valid value."))
            collection = choose_msigdb_collection(version)
        }
    }
    
    url = paste0(msigdb_base_url, "/", version, "/", collection, ".v", version, ".", gene_id_type, ".gmt")
    basename = basename(url)
    if(is.null(msigdb_env$gene_sets[[basename]])) {
        con = url(url)
        on.exit(close(con))
        ln = readLines(con)
        ln = strsplit(ln, "\t")
        gs = lapply(ln, function(x) x[-(1:2)])
        names(gs) = sapply(ln, function(x) x[1])
        msigdb_env$gene_sets[[basename]] = gs
    } else {
        gs = msigdb_env$gene_sets[[basename]]
    }
   
    if(as_table) {
        df = data.frame(gene_set = rep(names(gs), times = sapply(gs, length)),
                        gene = unlist(gs))
        rownames(df) = NULL
        df
    } else {
        gs
    }
}


map_to_entrez_id = function (from, org_db = "org.Hs.eg.db") {
    
    if(inherits(org_db, "character")) {
        org_db = getFromNamespace(org_db, ns = org_db)
    }

    x = keys(org_db, keytype = from)
    if("GENETYPE" %in% columns(org_db)) {
        map_tb = AnnotationDbi::select(org_db, keys = x, keytype = from, columns = c("ENTREZID", "GENETYPE"))
        map_tb = map_tb[map_tb$GENETYPE == "protein-coding", , drop = FALSE]
    } else {
        map_tb = AnnotationDbi::select(org_db, keys = x, keytype = from, columns = c("ENTREZID"))
    }
    map_tb = map_tb[!duplicated(map_tb[, 1]), , drop = FALSE]
    
    structure(map_tb[, 2], names = map_tb[, 1])
}

#' @importFrom utils getFromNamespace 
guess_id_type = function (id, org_db = "org.Hs.eg.db", verbose = TRUE) {
    l = grepl("^\\d+$", id)
    if (sum(l)/length(l) > 0.5) {
        return("ENTREZID")
    }
    l = grepl("^ENS.*G", id)
    if (sum(l)/length(l) > 0.5) {
        return("ENSEMBL")
    }
    l = grepl("^ENS.*T", id)
    if (sum(l)/length(l) > 0.5) {
        return("ENSEMBLTRANS")
    }
    l = grepl("^(NC|NG|NM|NR|NP|XM|XR|XP|WP)_\\d+", id)
    if (sum(l)/length(l) > 0.5) {
        return("REFSEQ")
    }

    if(inherits(org_db, "OrgDb")) {
        return(NULL)
    }

    if(inherits(org_db, "character")) {
        org_db = getFromNamespace(org_db, ns = org_db)
    }

    all_ids = keys(org_db, keytype = "SYMBOL")
    
    l = sample(id, min(100, length(id))) %in% all_ids
    p_match = sum(l)/length(l)
    if (p_match > 0.5) {
        if (verbose) cat("  gene id might be SYMBOL (p = ", sprintf('%.3f', p_match), ")\n")
        return("SYMBOL")
    }
    
    if (verbose) ("  cannot decide which gene id to use.\n")
    return(NULL)
}

guess_id_mapping = function (id, org_db = "org.Hs.eg.db", verbose = TRUE) {
    col = guess_id_type(id, org_db, verbose = verbose)
    if (is.null(col)) {
        return(NULL)
    }
    if (col == "ENTREZID") {
        return(NULL)
    }
    id_mapping = map_to_entrez_id(col, org_db)
    l = grepl("^ENS.*(G|T)", id) | grepl("^(NC|NG|NM|NR|NP|XM|XR|XP|WP)_\\d+", id)
    if (sum(l)/length(l) > 0.5) {
        fun = local({
            id_mapping = id_mapping
            function(x) {
                x = gsub("\\.\\d+$", "", x)
                id_mapping[x]
            }
        })
        return(fun)
    }
    else {
        return(id_mapping)
    }
}

#' Convert to EntreZ IDs
#' 
#' @rdname convert_to_entrez_id
#' 
#' @param x A vector of gene IDs, a named numeric vector where gene IDs are names, or a matrix with gene IDs as row names
#' @param org_db The name of the `OrgDb` object.
#' 
#' @export
convert_to_entrez_id = function(x, org_db = "org.Hs.eg.db") {
    if(is.matrix(x)) {
        map = guess_id_mapping(rownames(x), org_db = org_db)
        if(is.null(map)) stop("Cannot detect gene ID type.")
        if(is.function(map)) {
            new_rn = map(rownames(x))
        } else {
            new_rn = map[rownames(x)]
        }
        l = is.na(new_rn)

        x = x[!l, , drop = FALSE]
        new_rn = new_rn[!l]

        x2 = do.call(rbind, tapply(1:nrow(x), new_rn, function(ind) {
            colMeans(x[ind, , drop = FALSE])
        }))
        return(x2)

    } else if(is.numeric(x)) {
        map = guess_id_mapping(names(x), org_db = org_db)
        if(is.null(map)) stop("Cannot detect gene ID type.")
        x2 = x
        if(is.function(map)) {
            names(x2) = map(names(x))
        } else {
            names(x2) = map[names(x)]
        }
        x2 = x2[!is.na(names(x2))]
        x2 = tapply(x2, names(x2), mean)
        return(structure(as.vector(x2), names = names(x2)))
    } else {
        map = guess_id_mapping(x, org_db = org_db)
        if(is.null(map)) stop("Cannot detect gene ID type.")
        if(is.function(map)) {
            x2 = map(x)
        } else {
            x2 = map[x]
        }
        x2 = x2[!is.na(x2)]
        x2 = x2[!duplicated(x2)]
        return(x2)
    }
}


#' Get GO gene sets from OrgDb object
#' 
#' @param orgdb An `OrgDb` object.
#' @param ontology "BP", "CC" or "MF".
#' @param gene_id_type Which gene ID type to use.
#' @param as_table Whether to return a list or a two-column data frame?
#' 
#' @export
get_GO_gene_sets_from_orgdb = function(orgdb, ontology = "BP", 
    gene_id_type = c("ENTREZID", "SYMBOL", "ENSEMBL"), as_table = FALSE) {
    
    ontology = toupper(ontology)

    gene_id_type = match.arg(gene_id_type)[1]

    tb = AnnotationDbi::select(orgdb, keys = keys(orgdb, gene_id_type), columns = c("GOALL", "ONTOLOGYALL"), keytype = gene_id_type)
    tb = tb[tb$ONTOLOGYALL %in% ontology, , drop = FALSE]
    tb = tb[!is.na(tb$GOALL), c(gene_id_type, "GOALL"), drop = FALSE]
    tb = unique(tb)
    if(as_table) {
        tb
    } else {
        split(tb[, gene_id_type], tb$GOALL)
    }
}

#' Random genes
#' 
#' @param orgdb An `OrgDb` object.
#' @param n Number of random genes.
#' @param keytype Keytype of the genes in the `OrgDb` database.
#' 
#' @export
#' @import AnnotationDbi
random_genes = function(orgdb, n = 1000, keytype = "ENTREZID") {
    if("GENETYPE" %in% columns(orgdb)) {
        gi = AnnotationDbi::select(orgdb, key = "protein-coding", keytype = "GENETYPE", column = keytype)[, 2]
    } else {
        gi = keys(orgdb, keytype = keytype)
    }
    sample(gi, min(n, length(gi)))
}


#' Over-representation analysis
#' 
#' @param genes A vector of genes.
#' @param gene_sets A list of vectors.
#' @param universe A vector of background genes.
#' @param min_hits Minimal number of genes in `genes` in a gene set.
#' 
#' @importFrom stats phyper
#' @importFrom stats p.adjust
#' @export
ora = function(genes, gene_sets, universe = NULL, min_hits = 3) {

    if(is.null(universe)) {
        universe = unique(unlist(gene_sets))
    } else {
        universe = unique(universe)
    }
    # restrict in universe
    genes = intersect(genes, universe)
    gene_sets = lapply(gene_sets, function(x) intersect(x, universe))

    n_universe = length(universe)
    n_genes = length(genes)
    
    x = sapply(gene_sets, function(x) length(intersect(x, genes)))
    m = sapply(gene_sets, length)
    n = n_universe - m
    k = n_genes
    
    p = phyper(x - 1, m, n, k, lower.tail = FALSE)

    df = data.frame(
        gene_set = names(gene_sets),
        hits = x,
        n_genes = k,
        n_gs = m,
        n_total = n_universe,
        log2fe = log2(x*n_universe/k/m),
        p_value = p
    )
    df = df[df$hits >= min_hits, , drop = FALSE]
    df$p_adjust = p.adjust(df$p_value, "BH")
    rownames(df) = df$gene_set
    df[order(df$p_adjust, df$p_value), ,drop = FALSE]
}


#' Interactive Table for AnnotationHub
#' 
#' @param ah A `AnnotationHub` object.
#' 
#' @importFrom shiny fluidPage shinyApp
#' @importFrom htmltools h2 hr
#' @import DT
#' @import DBI
#' @import RSQLite
#' @export
ah_shiny = function(ah) {

    sqlite = dbDriver("SQLite")
    db = dbConnect(sqlite, ah@.db_path)

    tb = dbGetQuery(db, "SELECT r.ah_id, r.title, r.species, r.taxonomyid, r.genome, r.description, p.rdataclass, r.rdatadateadded from resources r inner join rdatapaths p on r.id = p.id")

    dbDisconnect(db)

    ui = fluidPage(
        title = 'AnnotationHub',
        h2("Interactive Table for AnnotationHub"),
        hr(),
        DT::dataTableOutput('tbl')
    )

    server = function(input, output, session) {
        output$tbl = DT::renderDataTable(tb, server = TRUE, filter = "top")
    }

    print(shinyApp(ui, server))
}


#' Add more columns
#' 
#' @param res A `enrichResult` object or a `gseaResult` object
#' 
#' @details
#' It adds the following columns in the enrichment table:
#' 
#' - n_hits
#' - n_genes
#' - gs_size
#' - n_total
#' - log2_fold_enrichment
#' - z_score
#' @export
add_more_columns = function(res) {
    k = as.numeric(gsub("/\\d+$", "", res@result$GeneRatio))
    m1 = as.numeric(gsub("^\\d+/", "", res@result$GeneRatio))
    m2 = as.numeric(gsub("/\\d+$", "", res@result$BgRatio))
    n = as.numeric(gsub("^\\d+/", "", res@result$BgRatio))

    res@result$n_hits = k
    res@result$n_genes = m1
    res@result$gs_size = m2
    res@result$n_totle = n

    res@result$log2_fold_enrichment = log2(k*n/m1/m2)
    res@result$z_score = (k - m1*m2/n)/sqrt(m1*m2*(n-m2)*(n-m1)/n/n/(n-1))

    res
}

#' Wrap long text into several lines
#' 
#' @param x A vector of sentences.
#' @param width Maximal number of chararacters per line.
#' @export
wrap_text = function(x, width = 60) {
    sapply(x, function(txt) {
        paste(strwrap(txt, width = width), collapse = "\n")
    })
}


#' GSEA analysis
#' 
#' @param expr The complete expression matrix
#' @param condition The condition labels of samples
#' @param cmp A vector of two, `cmp[1] - cmp[2] > 0` means up-regulation
#' @param geneset A vector of genes
#' @param plot Whether to make the GSEA plot?
#' @param power Power added to the weight.
#' 
#' @rdname GSEA
#' @export
#' @importFrom graphics abline points
#' @importFrom matrixStats rowSds
gsea = function(expr, condition, cmp, geneset, plot = FALSE, power = 1) {

    m1 = expr[, condition == cmp[1]]
    m2 = expr[, condition == cmp[2]]

    s = (rowMeans(m1) - rowMeans(m2))/(rowSds(m1) + rowSds(m2))

    s = sort(s, decreasing = TRUE)

    l_set = names(s) %in% geneset
    
    # f1 = cumsum(l_set)/sum(l_set)  # <<-- the original line
    
    s_set = abs(s)^power   # <<-- here
    s_set[!l_set] = 0
    f1 = cumsum(s_set)/sum(s_set)  ## <<- here

    l_other = !l_set
    f2 = cumsum(l_other)/sum(l_other)

    if(plot) {
        plot(f1 - f2, type = "l", xlab = "sorted genes")
        abline(h = 0, lty = 2, col = "grey")
        points(which(l_set), rep(0, sum(l_set)), pch = "|", col = "red")
        abline(v = which.max(f1 - f2), lty = 3, col = "blue")
    }

    max(f1 - f2)
}

#' @param s A pre-sorted gene-level scores.
#' @param perm Whether to perform gene-permutation?
#' 
#' @rdname GSEA
#' @export
gsea_gene_perm = function(s, geneset, perm = FALSE, power = 1) {
    
    if(perm) {
        # s is still sorted, but the gene labels are randomly shuffled
        # to break the associations between gene scores and gene labels.
        names(s) = sample(names(s))  ## <<- here
    }

    l_set = names(s) %in% geneset
    s_set = abs(s)^power
    s_set[!l_set] = 0
    f1 = cumsum(s_set)/sum(s_set)

    l_other = !l_set
    f2 = cumsum(l_other)/sum(l_other)

    max(f1 - f2)
}
