library(ROCR)
set.seed(1)

countdata = read.csv("/Users/elikond/Downloads/counts.cts", sep = '\t')
temp_col_names <- sub('.', '', colnames(countdata)[-c(1)])
new_col_names <- c("Genes",temp_col_names)
colnames(countdata) <- new_col_names
countdata <- countdata[!duplicated(countdata$Genes), ]
countdata <- na.omit(countdata)
rownames(countdata) <- countdata$Genes
labelled_countdata <- countdata[-c(1,10,12,14,15,17,21,27,28)]

coldata = read.csv("/Users/elikond/Downloads/Engraftment - 1.csv", sep = ",")
#rownames(coldata) <- sub('-', '.', coldata$Patient)
#coldata <- coldata[-c(1)]

keep <- which(coldata$Engraftment == 'h')
#Set 1 = High
high_labelled_counts <- labelled_countdata[,keep]
#Set 2 = Low
low_labelled_counts <- labelled_countdata[,-c(keep)]

#Finding column with least variance to remove so number of labelled and unlabelled patients is equal
var_low_labelled_counts <- colVars(as.matrix(low_labelled_counts[sapply(low_labelled_counts, is.numeric)]))
low_var_indx <- which(var_low_labelled_counts == min(var_low_labelled_counts))
low_labelled_counts <- low_labelled_counts[-c(low_var_indx)]

#setequal(rownames(low_labelled_counts),rownames(high_labelled_counts))
#[1] TRUE

mahalanobis_low_high <- function(high_distribution, low_distribution){
  distribution_list <- list(high_distribution, low_distribution)
  i = 0
  for (distribution in distribution_list) {
    i <- i + 1
    #Distance from high distribution
    distance_high = mahalanobis(distribution, colMeans(high_distribution), cov(high_distribution))
    pval_high = pchisq(distance_high, df=3, lower.tail=FALSE)
    
    #Distance from low Distribution
    distance_low = mahalanobis(distribution, colMeans(low_distribution), cov(low_distribution))
    pval_low = pchisq(distance_low, df=3, lower.tail=FALSE)
    
    if (i == 1){
      type <- rep(c(1),times=nrow(high_distribution))
      final_high_df <- data.frame(distance_high, pval_high, distance_low, pval_low, type)
    } else{
      type <- rep(c(0),times=nrow(low_distribution))
      final_low_df <- data.frame(distance_high, pval_high, distance_low, pval_low, type)
    }
  }
  final_df <- rbind(final_high_df, final_low_df)
  return(final_df)
}

df <- mahalanobis_low_high(high_labelled_counts, low_labelled_counts)
df$type <- as.integer(as.logical(df$type))

#Splitting df into train and test
split1<- sample(c(rep(0, 0.7 * nrow(df)), rep(1, 0.3 * nrow(df))))
train <- na.omit(df[split1 == 0, ])
test <- na.omit(df[split1== 1, ])

#Running model and summary stats
model <- glm(type ~.,family=binomial(link='logit'), data=train)
summary(model)

#Finding accuracy
fitted.results <- predict(model,newdata= test,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test$type)
print(paste('Accuracy',1-misClasificError))

#Finding AUC
p <- predict(model, newdata=subset(test), type="response")
pr <- prediction(p, test$type)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
