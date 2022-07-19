dirn <- system.file('extdata', 'GeneCalling', 'IdTaxaTraining',
                    package='CompGenomicsBioc2022')


trainingSet <- vector('list', 16)
for ( i in 1:16 ){
  f <- file.path(dirn, paste0('Training',i,'.RData'))
  if (i!=9 && file.exists(f)){
    load(f)
    trainingSet[[i]] <- dataentry
  }
}

ninthentry <- list()
for ( i in 1:10 ){
  f <- file.path(dirn, paste0('Training9_', i, '.RData'))
  load(f)
  ninthentry <- append(ninthentry, dataentry)
}
trainingSet[[9]] <- ninthentry

load(file.path(dirn, 'TrainingKEY.RData'))
names(trainingSet) <- trainingSetNames

class(trainingSet) <- c('Taxa', 'Train')

dirout <- system.file('extdata', 'GeneCalling',
                    package='CompGenomicsBioc2022')
outfile <- file.path(dirout, 'IdTaxaActinobacteriaTrainingSet.RData')
save(trainingSet, file=outfile)