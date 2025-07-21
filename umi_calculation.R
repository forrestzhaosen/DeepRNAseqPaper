
calculate_umi <- function(depth,umi){
  not_selected <- ((4^umi-1)/4^umi)^depth
  select_once <- depth*((4^umi-1)/4^umi)^(depth-1)/4^umi
  return((1-not_selected-select_once)*4^umi)
}

for (length in 5:12){
  calculate_umi(1000,length)
}

for (depth in c(10000,100000,1000000)){
  print(calculate_umi(depth,10)/depth)
}

X=5000000000
Y=3000000000

1-X*(1-2.71828^(-Y/X))/Y
