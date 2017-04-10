namedList = function(elements, names) {
  l = list(elements)
  names(l) = names
  l
}

getSubsets = function(elements, proper = TRUE) {
  ret = lapply(0:(length(elements) - ifelse(proper, 1, 0)),
    function(m) combn(elements, m))
  unlist(ret, FALSE)
}