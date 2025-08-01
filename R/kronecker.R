age <- matrix(1:100, 10, 10)
social_activity <- matrix((1:4) / 4, 2, 2)
age_social_activity <- kronecker(social_activity,
                                 age,
                                 FUN = "*")
par(mfrow = c(1, 3))
col_res <- 100
image(age,
      col = hcl.colors(col_res, "YlOrRd", rev = TRUE),
      asp = 1,
      axes = FALSE,
      main = "age")
image(social_activity,
      col = hcl.colors(col_res, "YlOrRd", rev = TRUE),
      asp = 1,
      axes = FALSE,
      main = "social activity")
image(age_social_activity,
      col = hcl.colors(col_res, "YlOrRd", rev = TRUE),
      asp = 1,
      axes = FALSE,
      main = "age * social activity")
