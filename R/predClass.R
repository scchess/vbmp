`predClass` <-
function(obj) {
    as.numeric(apply(obj$Ptest, 1, which.max));
}

