# install.packages("reader") if not extant
# Clear variables
rm(list=ls())

# Import reader to read select lines
library("reader")

############# Getting a, b, c, alpha, beta, gamma ##########################

# Reading the lines containing info from cif file (lines 68-73)
raw_data = n.readLines(paste("293146.cif", sep = ""),
            header = FALSE, skip = 67, n = 6)
# Converting the lines into a table
raw_data_table = read.table(textConnection(raw_data))
# Removing everything in brackets
raw_data_table$V2 = gsub("\\s*\\([^\\)]+\\)","",as.character(raw_data_table$V2))
# Extracting only numbers (from column 2) as a data frame
data_extracted = as.data.frame(raw_data_table[,2], drop=FALSE)
# Converting the data frame into matrix
data = as.matrix(data_extracted)
# Converting the matrix into a numeric vector
data = as.numeric(data)

a = data[1]
b = data[2]
c = data[3]
alpha = data[4]
beta = data[5]
gamma = data[6]

################ Getting fractional coordinates ######################
# Reading the lines containing info from cif file (lines 187-225)
raw_coord = n.readLines(paste("293146.cif", sep = ""),
                       header = FALSE, skip = 186, n = 39)
# Converting the lines into a table
raw_coord_table = read.table(textConnection(raw_coord))
# Removing everything in brackets
raw_coord_table$V3 = gsub("\\s*\\([^\\)]+\\)","",as.character(raw_coord_table$V3))
raw_coord_table$V4 = gsub("\\s*\\([^\\)]+\\)","",as.character(raw_coord_table$V4))
raw_coord_table$V5 = gsub("\\s*\\([^\\)]+\\)","",as.character(raw_coord_table$V5))
# Extracting only numbers (from columns 3 to 5) as a data frame
x_coord_extracted = as.data.frame(raw_coord_table[,3], drop=FALSE)
y_coord_extracted = as.data.frame(raw_coord_table[,4], drop=FALSE)
z_coord_extracted = as.data.frame(raw_coord_table[,5], drop=FALSE)
atoms = as.data.frame(raw_coord_table[,2], drop=FALSE)
# Converting the data frame into matrix
x_coord = as.matrix(x_coord_extracted)
y_coord = as.matrix(y_coord_extracted)
z_coord = as.matrix(z_coord_extracted)
atoms = as.matrix(atoms)
colnames(atoms) <- "Atoms"
# Converting the matrix into a numeric vector
x_coord = as.numeric(x_coord)
y_coord = as.numeric(y_coord)
z_coord = as.numeric(z_coord)

coord = rbind(x_coord, y_coord, z_coord)

######### Converting fractional to orthogonal ################
cos_alpha = round(cos(alpha*(pi/180)), digits = 4)
cos_beta = round(cos(beta*(pi/180)), digits = 4)
cos_gamma = round(cos(gamma*(pi/180)), digits = 4)
sin_gamma = round(sin(gamma*(pi/180)), digits = 4)

v = round(sqrt(1 - (cos_alpha^2) - (cos_beta^2) - (cos_gamma^2) + (2 * cos_alpha * cos_beta * cos_gamma)), digits = 4)
conversion_matrix = matrix(c(a, (b * cos_gamma), (c * cos_beta), 0, (b * sin_gamma), (c * (cos_alpha - (cos_beta * cos_gamma)) / sin_gamma),
                             0, 0, (c * v / sin_gamma)),
                           nrow = 3, ncol = 3, byrow = TRUE)
orth = matrix(rep(0, each = 3), nrow = ncol(coord), ncol = 3)
colnames(orth) <- c("x", "y", "z")
for(i in 1:ncol(coord)){
  orth[i,] = round(conversion_matrix %*% coord[,i], digits = 3)
}
pdb = cbind.data.frame(atoms,orth)
# write.table(pdb, file = "mypdb.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
