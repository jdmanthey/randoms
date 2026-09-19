options(scipen=999)
require(GenomicRanges)

input_file <- "GCA_018697195.1_ASM1869719v1_genomic.fna.out"
output_file <- "GCA_018697195.1_ASM1869719v1_genomic.fna.no_dups.out"

# read repeatmasker file
x <- read.table(input_file, sep="", stringsAsFactors=F, fill=T, skip=2, col.names=c("SW_score", "perc_div", "perc_del", "perc_ins", "query", "q_start", "q_end", "q_leftover", "orientation", "repeat_match", "class/family", "r_start", "r_end", "r_leftover", "ID", "overlap"))

# all query sequences in genome (scaffolds)
x_chr <- unique(x$query)

# start output
write(colnames(x), file=output_file, ncolumns=16, sep="\t")

for(b in 1:length(x_chr)) {
	b_rep <- x[x$query == x_chr[b],]
	
	# find overlaps
	x_temp <- data.frame(chr=as.character(b_rep[,5]), start=as.numeric(b_rep[,6]), end=as.numeric(b_rep[,7]))
	overlaps <- findOverlaps(makeGRangesFromDataFrame(x_temp), makeGRangesFromDataFrame(x_temp))
	# convert to matrix and remove self matches
	overlaps <- cbind(overlaps@from, overlaps@to)
	overlaps <- overlaps[overlaps[,1] != overlaps[,2], ]
	# which rows are involved
	TEs_to_check <- sort(unique(c(overlaps[,1], overlaps[,2])))
	
	# populate bp for each TE to check
	bp_list <- list()
	for(d in TEs_to_check) {
		bp_list[[d]] <- seq(from=b_rep$q_start[d], to=b_rep$q_end[d])
	}
	# go through comparisons for each TE
	for(d in TEs_to_check) {
		d_overlaps <- overlaps[overlaps[,1] == d | overlaps[,2] == d,]
		d_overlaps <- sort(unique(c(d_overlaps[,1], d_overlaps[,2])))
		d_overlaps <- d_overlaps[d_overlaps != d]
		for(e in 1:length(d_overlaps)) {
			if(b_rep$perc_div[d] > b_rep$perc_div[d_overlaps[e]]) { # rm bp from list if more divergent than comp
				# remove bp found in other element
				bp_list[[d]] <- bp_list[[d]][bp_list[[d]] %in% bp_list[[d_overlaps[[e]]]] == F]
			} else if (b_rep$perc_div[d] == b_rep$perc_div[d_overlaps[e]]) { # check if equal div
				if(b_rep$SW_score[d] < b_rep$SW_score[d_overlaps[e]]) { # check if worse SW score
					# remove bp found in other element
					bp_list[[d]] <- bp_list[[d]][bp_list[[d]] %in% bp_list[[d_overlaps[[e]]]] == F]
				} else if(b_rep$SW_score[d] == b_rep$SW_score[d_overlaps[e]]) { # if equal SW score
					if(d > d_overlaps[e]) {
						# remove bp found in other element
						bp_list[[d]] <- bp_list[[d]][bp_list[[d]] %in% bp_list[[d_overlaps[[e]]]] == F]
					} else {
						# do nothing if this is first element and has same SW score
					}
				} else {
					# do nothing if this is equal divergent but has better SW score
				}
			} else {
				# do nothing if this is less divergent or has better SW score
			}
		}
	}
	
	# replace beginning and end in b_rep for each TE 
	extra_windows <- c()
	mark_for_removal <- c()
	for(d in TEs_to_check) {
		d_rep <- bp_list[[d]]
		if(length(d_rep) > 1) { # at least two bp in TE
			# check that all bp are 1 bp apart
			if(max(diff(d_rep)) == 1) {
				# replace start and end values of TE (may or may not change original value)
				b_rep$q_start[d] <- min(d_rep)
				b_rep$q_end[d] <- max(d_rep)
			} else { # if the bp are split apart
				d_diff <- diff(d_rep)
				# check that the split is not just the first or last bp
				if(d_diff[1] > 1) {
					d_rep <- d_rep[2:length(d_rep)]
				} else if(d_diff[length(d_diff)] > 1) {
					d_rep <- d_rep[1:(length(d_rep) - 1)]
				}
				d_diff <- diff(d_rep)
				if(max(d_diff) == 1) {
					# replace start and end values of TE (may or may not change original value)
					b_rep$q_start[d] <- min(d_rep)
					b_rep$q_end[d] <- max(d_rep)
				} else {
					d_diff <- c(1, d_diff)
					n_loops <- length(d_diff[d_diff > 1]) + 1
					d_breaks <- seq(from=1, to=length(d_diff))[d_diff > 1]
					# loop for each section
					for(e in 1:n_loops) {
						if(e == 1) {
							e_start <- d_rep[1]
							e_end <- d_rep[d_breaks[e] - 1]
							# replace values in b_rep
							b_rep$q_start[d] <- e_start
							b_rep$q_end[d] <- e_end
						} else if(e == n_loops){
							e_start <- d_rep[d_breaks[e-1]]
							e_end <- d_rep[length(d_rep)]
							# add to extra windows to be added
							row_template <- b_rep[d,]
							row_template$q_start <- e_start
							row_template$q_end <- e_end
							extra_windows <- rbind(extra_windows, row_template)
						} else {
							e_start <- d_rep[d_breaks[e-1]]
							e_end <- d_rep[d_breaks[e] - 1]
							# add to extra windows to be added
							row_template <- b_rep[d,]
							row_template$q_start <- e_start
							row_template$q_end <- e_end
							extra_windows <- rbind(extra_windows, row_template)
						}
					}
				}
			}
		} else {
			# if no bp or only a single bp, mark for removal
			mark_for_removal <- c(mark_for_removal, d)
		}	
	}
	
	# remove lines marked for removal
	if(length(mark_for_removal) > 0) { b_rep <- b_rep[-mark_for_removal,] }
	
	# add extra windows
	if(is.null(extra_windows) == FALSE) { b_rep <- rbind(b_rep, extra_windows) }
	
	# sort
	if(is.null(extra_windows) == FALSE) { b_rep <- b_rep[order(b_rep$q_start),] }

	# output to file
	write.table(b_rep, file=output_file, sep="\t", quote=F, row.names=F, col.names=F, append=T)
	
	
	# find overlaps
	x_temp <- data.frame(chr=as.character(b_rep[,5]), start=as.numeric(b_rep[,6]), end=as.numeric(b_rep[,7]))
	overlaps <- findOverlaps(makeGRangesFromDataFrame(x_temp), makeGRangesFromDataFrame(x_temp))
	# convert to matrix and remove self matches
	overlaps <- cbind(overlaps@from, overlaps@to)
	overlaps <- overlaps[overlaps[,1] != overlaps[,2], ]
	if(nrow(overlaps) > 0) {
		stop()
	}
}	
		
	
	
	
