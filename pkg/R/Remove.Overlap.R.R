`Remove.Overlap.R` <-
function(chpts){

    lookup <- nrow(chpts)
    
    ret.v <- rep(FALSE, lookup)
    ret.v[1] <- TRUE
    for(i in 2:lookup){
        overlap <- 0;
        st <- chpts[i,1];
        ed <- chpts[i,2];
        for(j in which(ret.v)){
            # check for overlap.
            if(     (st<=chpts[j,1] && ed>=chpts[j,2]) || 
                (st>=chpts[j,1] && st<=chpts[j,2]) || 
                (ed>=chpts[j,1] && ed<=chpts[j,2]) )
            {
                overlap <- 1;
                break;
            }
        }
        if(overlap == 0){
            ret.v[i] <- TRUE
        }
            
    }

    return(ret.v)
}

