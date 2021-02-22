s2n <- function(row, labels){
  control_row = row[labels$Group=="c"]
  disease_row = row[labels$Group=="d"]
  ratio <- (mean(disease_row)-mean(control_row))/(sd(control_row)+sd(disease_row))
  ratio
}


lfc <- function(row, labels){
  control_row = row[labels$Group=="c"]
  disease_row = row[labels$Group=="d"]
  ratio <- log(mean(disease_row)/mean(control_row))
  ratio
}

rank_genes <- function(data, labels, rank='s2n'){
  if (rank == 's2n'){
    ranks <- apply(data, 1, s2n, labels)
  }else if(rank=='lfc'){
    ranks <- apply(data, 1, lfc, labels)
  }
  ranks
}
