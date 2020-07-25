# Load the library multtest, ROCR and Load golub data 
library(multtest)
library(ROCR)
data(golub)

# Question i. Get the Size of the dataset
size.of.the.dataset = dim(golub)
cat("Question i. Size of the dataset: ", size.of.the.dataset, "\n")
cat("\n")

# Question ii. View one of the gene names - use row and column numbers. My name as the variable name
geneNames = golub.gnames[1042,2]
cat("Question ii. Name of the gene 1042: ", geneNames, "\n")
cat("\n")

# Question iii. View the expression profiles of this gene across all patients
expression.profiles = golub[1042,]
cat("Question iii. The expression profiles of gene 1042 across all patients: ", expression.profiles, "\n")
cat("\n")

# Question iv. View the class labels
class.labels = golub.cl
cat("Question iv. The class labels: ", class.labels, "\n")
cat("\n")

# Question v. How many patients are there in each class
no.of.patients.in.class.zero = length(which(golub.cl == 0)) # Number of Patients in Class Zero
no.of.patients.in.class.one = length(which(golub.cl == 1)) # Number of Patients in Class One
cat(sprintf("Question v. Number of Patients in Class Zero = %i \n\t\t Number of Patients in Class One = %i \n", no.of.patients.in.class.zero, no.of.patients.in.class.one))
cat("\n")

# Question vi. Change labels to ALL and AML (0/1 not convenient)
new.class.labels = factor(golub.cl, levels = c(0,1), labels = c("ALL","AML"))
cat("Question vi. New class labels: ", "\n")
print(new.class.labels)
cat("\n")

# Question vii. View expressions of gene 1042 only for ALL class
expressions.for.ALL.class = golub[1042, new.class.labels == "ALL"]
cat("Question vii. Expressions of gene 1042 only for ALL class: ", expressions.for.ALL.class, "\n")
cat("\n")

# Question viii. View expressions of gene 1042 only for AML class
expressions.for.AML.class = golub[1042, new.class.labels == "AML"]
cat("Question vii. Expressions of gene 1042 only for AML class: ", expressions.for.AML.class, "\n")
cat("\n")

# Question ix. View the expression values for the two groups of this gene as boxplots. Is this a good discriminant gene?
boxplot(golub[1042, new.class.labels == "ALL"], golub[1042, new.class.labels == "AML"], names = c("ALL","AML"), xlab = "Class", ylab = "Expression Values")
cat("Question ix. The expression values for the two groups of gene 1042 as boxplots: ", "\n")
cat("\n")
cat("\t", "Is this a good discriminant gene? ", "According to the boxplot results for the two groups of gene 1042, the two groups have very little overlap. Thus, gene 1042 can be considered as a good discriminant gene.", "\n")
cat("\n")

# Question x. Try a few random genes and indicate if they are discriminant or not.  
cat("Question x. Try a few random genes and indicate if they are discriminant or not: Tested Genes: 50, 200, 267, 600, 2100, 2420, 2800, 3000, 3020, 3042", "\n")
cat("\t", "According to the tested random genes, most of the genes overlap too much and do not discriminate the two groups well compared to gene 1042.", "\n")
cat("\n")

# Question xi. Compute the mean of ALL subjects for gene 1042.  My name as the variable name.
MeanALLSubjects = mean(golub[1042, new.class.labels == "ALL"])
cat("Question xi. The mean of ALL subjects for gene 1042: ", MeanALLSubjects, "\n")
cat("\n")

# Question xii. Compute the mean of AML for the same gene.
MeanAML = mean(golub[1042, new.class.labels == "AML"])
cat("Question xii. The mean of AML for gene 1042: ", MeanAML, "\n")
cat("\n")

# Question xiii. Compute the Fisher Ratio for this gene
VarienceALL = var(golub[1042, new.class.labels == "ALL"])
VarienceAML = var(golub[1042, new.class.labels == "AML"])
fisher.ratio = (MeanALLSubjects - MeanAML)^2 / (VarienceALL + VarienceAML)
cat(sprintf("Question xiii. Fisher Ratio for gene 1042: %.4f\n", fisher.ratio))
cat("\n")

# Question xiv. Compute Fisher Ratios of all the genes. 
no.of.genes = nrow(golub)
fisher.ratio.vector = vector()
cat("Question xiv. Fisher Ratios of all the genes: ", "\n")
for (gene.id in 1:no.of.genes) {
  current.MeanALLSubjects = mean(golub[gene.id, new.class.labels == "ALL"])
  current.MeanAML = mean(golub[gene.id, new.class.labels == "AML"])
  current.VarienceALL = var(golub[gene.id, new.class.labels == "ALL"])
  current.VarienceAML = var(golub[gene.id, new.class.labels == "AML"])
  current.fisher.ratio = (current.MeanALLSubjects - current.MeanAML)^2 / (current.VarienceALL + current.VarienceAML)
  cat("\t", sprintf("Fisher Ratio for gene %i: %f\n", gene.id, current.fisher.ratio))
  fisher.ratio.vector[gene.id] = current.fisher.ratio
}
cat("\n")

# Question xv. How are these Fisher ratios distributed? Show using a histogram.  Can you identify the discriminant genes?
cat("Question xv. Distribution of Fisher Ratios using a histogram: ", "\n")
hist(fisher.ratio.vector)
cat("\n")
cat("\t", "Can you identify the discriminant genes? ", "", "\n")
cat("\t","\t","In general good discriminant genes have higher Fisher Ratios. According to the above histogram, most of genes have low fisher ratios and there are only a few genes with higher fisher ratios. Accordingly, there is a low frequency of significantly discriminant genes in this dataset.", "\n")
cat("\n")

# Question xvi. The article by Golub says CD33 plays an important role  in discriminating these two types of leukhemia  Which is it?
cat("Question xvi. CD33 plays an important role  in discriminating these two types of leukhemia  Which is it? ", "Acute Myeloid Leukemia and Acute Lymphoblastic Leukemia", "\n")
cat("\n")

# Question xvii. Collect the 10 genes for which Fisher Ratio is largest
cat("Question xvii. The 10 genes for which Fisher Ratio is largest: ", "\n")
temp = order(fisher.ratio.vector, decreasing = T)[1:10]
for(i in 1:length(temp)) {
 	cat("\t", i, "\t", golub.gnames[temp[i], 2], "\n")
}
cat("\n")

# Question xviii.	Collect expression levels of 50 most discriminant genes and view them as an intensity plot.
cat("Question xviii. Expression levels of 50 most discriminant genes and view them as an intensity plot: ", "\n")
top.fifty.discriminant.genes = sort(fisher.ratio.vector, decreasing=TRUE)
top.fifty.discriminant.genes = top.fifty.discriminant.genes[1:50]
den = density(top.fifty.discriminant.genes)
plot(den)
cat("\n")

# Question xix. Hierarchical cluster golub patients
cat("Question xix. Hierarchical cluster golub patients: ", "\n")
hierarchical.cluster.golub.data <- data.frame(golub[,]);
plot(hclust(dist(hierarchical.cluster.golub.data)));
cat("\n")

# Question xx. Plot ROC curves for classifying on genes CCND3 Cyclin D3 and gene Gdf5
cat("Question xx. ROC curves for classifying on genes CCND3 Cyclin D3 and gene Gdf5: ", "\n")
true.label <- factor(golub.cl, levels = 0:1, labels = c("TRUE","FALSE"));
pred.one <- prediction(golub[which(golub.gnames[,2] == "CCND3 Cyclin D3"),], true.label);
pred.two <- prediction(golub[which(golub.gnames[,2] == "Gdf5 gene"),], true.label);
perf.one <- performance(pred.one, "tpr", "fpr");
perf.two <- performance(pred.two, "tpr", "fpr");
plot(perf.one, col = "tomato", main="ROC curves for classifying on genes CCND3 Cyclin D3 and gene Gdf5");
plot(perf.two, add = TRUE, col = "turquoise");
cat("\n")






