subjectNum = 12
chanIntVec = c(6)
chanIntConditioningPairVec = c(1)
blockIntLM = c(1,3,5,6,7)
blockIntPlot = c(1,3,5,6,7,8)
blockNames <- c(
  `1` = "baseline 1",
  `3` = "baseline 2",
  `5` = "post A/B 200 ms 1",
  `6` = "baseline - post A/B 200 ms 2/pre A/A",
  `7` = "post A/A 200 ms 1",
  `8` = "baseline - post A/A 200 ms 2"
)
blockType = c('baseline','baseline','A/B 200',
              'baseline','A/A 200','A/A 200')
#whichCompareVec = list(c(1,3),c(3,5),c(6,7))
whichCompareVec = list(c(3,5),c(6,7))
