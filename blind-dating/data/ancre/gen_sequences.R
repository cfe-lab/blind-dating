library(RSQLite)
library(ape)

source('include/raxml.R')

sql     <- dbDriver("SQLite")
ancredb <- dbConnect(sql, "ancre.sqlite")

patients <- dbGetQuery(ancredb, "SELECT DISTINCT p.PatientID, p.TimeType, p.Treated, p.SuperInfection FROM Patients AS p 
	JOIN Sequences AS s ON p.PatientID == s.PatientID 
	JOIN Products AS pr ON pr.Accession == s.Accession
	JOIN SequenceType AS t ON t.Accession == s.Accession 
	WHERE t.isDNA == 0 AND 
		(p.TimeType == 'DPI' or p.TimeType == 'DPSC') AND 
		(pr.Product == 'env') AND p.SuperInfection == 0 AND p.Treated == 0;")


get.labels <- function(patient, is.dna){

	query <- sprintf( 
"SELECT DISTINCT pr.Accession, pr.Product, p.TimeType, s.days FROM Patients AS p 
	JOIN Sequences AS s ON p.PatientID == s.PatientID 
	JOIN Products AS pr ON pr.Accession == s.Accession
	JOIN SequenceType AS t ON t.Accession == s.Accession 
	WHERE  
		(p.TimeType == 'DPI' or p.TimeType == 'DPSC') AND 
		(pr.Product == 'env') AND p.SuperInfection == 0 AND 
		p.Treated == 0 AND
		t.isDNA == '%d' AND p.PatientID == '%d';", as.integer(is.dna), as.integer(patient)) 
	 dbGetQuery(ancredb, query)
}

for(p in patients$PatientID) { 
	early <- sprintf('1_align/%d_env_early.fasta', p)
	late <- sprintf('1_align/%d_env_late.fasta', p)

	rna_labels <- dbGetQuery(ancredb, "SELECT DISTINCT p.PatientID, p.TimeType, p.Treated, p.SuperInfection FROM Patients AS p 
	JOIN Sequences AS s ON p.PatientID == s.PatientID 
	JOIN Products AS pr ON pr.Accession == s.Accession
	JOIN SequenceType AS t ON t.Accession == s.Accession 
	WHERE t.isDNA == 0 AND 
		(p.TimeType == 'DPI' or p.TimeType == 'DPSC') AND 
		(pr.Product == 'env') AND p.SuperInfection == 0 AND p.Treated == 0;")
	# print(p) 

	rna.labs <- apply( get.labels(p, 0), 1, function(row) sprintf("%s_%s_%s_%d",row[1], row[2], row[3], as.integer(row[4])))
	dna.labs <- apply( get.labels(p, 1), 1, function(row) sprintf("%s_%s_%s_%d",row[1], row[2], row[3], as.integer(row[4])))

	# print(rna.labs)
	# print(dna.labs)

	early <- read.FASTA(early)
	late <- read.FASTA(late)

	dna <- rbind(rbind(as.matrix(late), as.matrix(early)))
	write.dna(dna, sprintf('fasta_wo_ref/patient_%d.fasta', p), format="fasta")
	selected <- c()

	if(length(rna.labs) && length(dna.labs)){
		print(sprintf("This patient %d has both!", p))
	}

	for(rna in rna.labs) {
		if(rna %in% labels(dna)) {
			seq <- as.matrix(dna[which(labels(dna)== rna), ])
			if(is.null(selected)){
				selected <- seq
			} else {
				selected <- rbind(selected, seq)
			}
		}
	}
#	write.dna(selected, sprintf('output_%d.whatever.seq', p))

	# tree <- raxml(selected, N=100, parsimony.seed=10000, bootstrap.seed=1000)
	# write.tree(tree, sprintf("7_ml_tree/%d_env_ml_rna.nwk", p))
}


