subjectNum = 20
chanIntVec = c(5)
chanIntConditioningPairVec = c(1)
blockIntLM = c(1,2,3,4,5,6,7,8,9)

blockIntPlot = blockIntLM

blockNames <- c(
  `1` = "baseline 1",
  `2` = "baseline 2 (pre conditioning)",
  `3` = "post A/B 200 ms (5 minutes)",
  `4` = "baseline 3 (post 200 ms A/B)",
  `5` = "baseline 4, after clinical testing",
  `6` = "post A/B 200 ms (15 mins)",
  `7` = "baseline 5 (post A/B)",
  `8` = "post A/A 200 ms (15 mins)",
  `9` = "baseline 6"
  )

blockType = c('baseline','baseline','A/B 200 5','baseline','baseline',
              'A/B 200','baseline','A/A 200','baseline')

whichCompareVec = list(c(5,6),c(7,8))

