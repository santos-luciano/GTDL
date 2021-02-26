## code to prepare `artset1987` dataset goes here

usethis::use_data(artset1987, overwrite = TRUE)

artset1987<-c(0.1,0.2,1,1,1,1,1,2,3,6,7,11,12,18,18,18,18,18,21,32,36,40,45,46,47,50,55,60,63,63,67,67,67,67,72,75,79,82,82,83,84,84,84,85,85,85,85,86,86)


write_csv(artset1987,"data-raw/artset1987.csv")
save(artset1987,file = "data/artset1987.rda")