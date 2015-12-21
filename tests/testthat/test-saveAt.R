data(wheat)
attach(wheat.Y)

if(.Platform$OS.type=="windows"){
  test_that("default save path",{lm1=FW(y,VAR,ENV,nIter=100,burnIn=50)})
  test_that("save path with no trailing file separation",{
  saveAt="C://tryFW"
    lm1=FW(y,VAR,ENV,nIter=100,burnIn=50,saveAt=saveAt);
    expect_equal(file.exists(paste(saveAt,"samps.rda",sep="")),T);
  }
  )
  test_that("save path with trailing file separation",{
  saveAt="C://tryFW//"
    lm2=FW(y,VAR,ENV,nIter=100,burnIn=50,saveAt=saveAt)
  expect_equal(file.exists(paste(saveAt,"samps.rda",sep="")),T)
  }
  )
}

#test save At on both unix and windows
test_that("default save path",{lm1=FW(y,VAR,ENV,nIter=100,burnIn=50)})
test_that("save path with no trailing file separation",{
  saveAt="~/tryFW"
  lm1=FW(y,VAR,ENV,nIter=100,burnIn=50,saveAt=saveAt);
  expect_equal(file.exists(paste(saveAt,"samps.rda",sep="")),T);
}
)
test_that("save path with trailing file separation",{
  saveAt="~/tryFW/"
  lm2=FW(y,VAR,ENV,nIter=100,burnIn=50,saveAt=saveAt)
expect_equal(file.exists(paste(saveAt,"samps.rda",sep="")),T)
}
)