ExtractSubSample <- function(data,boostrap_index,outofbag_index){
  bootstrap <- data
  outofbag <- data
  
  for(i in 1:length(data)){
    if(is.null(nrow(data[[i]]))){
      bootstrap_names <- names(data[[i]])[boostrap_index]
      bootstrap[[i]] <- data[[i]][boostrap_index]
      if(!is.null(bootstrap_names)){
        if(length(unique(bootstrap_names)) != length(bootstrap_names)){
          names(bootstrap[[i]]) <- make.names(bootstrap_names,unique = T)
        }
      }
      outofbag[[i]] <- data[[i]][outofbag_index]
    }else{
      sampleDim <- which.max(dim(data[[i]]))
      if(sampleDim == 1){
        bootstrap_names <- rownames(data[[i]])[boostrap_index]
        bootstrap[[i]] <- data[[i]][boostrap_index,]
        if(!is.null(bootstrap_names)){
          if(length(unique(bootstrap_names)) != length(bootstrap_names)){
            rownames(bootstrap[[i]]) <- make.names(bootstrap_names,unique = T)
          }
        }
        outofbag[[i]] <- data[[i]][outofbag_index,]
      }else{
        bootstrap_names <- colnames(data[[i]])[boostrap_index]
        bootstrap[[i]] <- data[[i]][,boostrap_index]
        if(!is.null(bootstrap_names)){
          if(length(unique(bootstrap_names)) != length(bootstrap_names)){
            colnames(bootstrap[[i]]) <- make.names(bootstrap_names,unique = T)
          }
        }
        outofbag[[i]] <- data[[i]][,outofbag_index]
      }
    }
  }
  return(list(bootstrap = bootstrap,outofbag = outofbag))
}
