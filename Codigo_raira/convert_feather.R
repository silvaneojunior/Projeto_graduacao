library(feather)

for(idade in c(1:5)){
for(i in c(1,6,12)){
  data=read.csv(paste0('Varicela/grid_data/MRE_pred_sem_delay__',i,'_',idade,'.csv'))
  data$Error=data$Error
  write_feather(data,paste0('Varicela/grid_data/MRE_pred_',i,'_',idade,'.feather'))
  write.csv(data,paste0('Varicela/grid_data/MRE_pred_',i,'_',idade,'.csv'))
}
}

data=read_feather(paste0('Varicela/grid_data/MRE_pred__',1,'.feather'))
