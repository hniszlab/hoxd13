
images=all_constructs_filtered_intensity
initial <- filter(images, images[,17] == 1) 

ml=list()
ml[[1]]=as.numeric(initial[which(initial[,19]=="wt"),5])
ml[[4]]=as.numeric(initial[which(initial[,19]=="DEdel"),5])
ml[[3]]=as.numeric(initial[which(initial[,19]=="-7A"),5])
ml[[2]]=as.numeric(initial[which(initial[,19]=="-15A"),5]) 

cor.test(ml[[1]],ml[[2]])
cor.test(ml[[1]],ml[[3]])
cor.test(ml[[1]],ml[[4]])
cor.test(ml[[2]],ml[[3]])
cor.test(ml[[2]],ml[[4]])
cor.test(ml[[3]],ml[[4]])
