subjectNum = 14
chanIntVec = c(5,8)
chanIntConditioningPairVec = c(1,0)
blockIntLM = c(1,2,3,4,5,6,7,8,9,10)

blockIntPlot = blockIntLM

blockNames <- c(
  `1` = "baseline 1",
  `2` = "baseline 2 (pre conditioning)",
  `3` = "post A/A 200 ms",
  `4` = "baseline 3 (post 200 ms)",
  `5` = "post A/B 200 ms",
  `6` = "baseline 4 (post 200 ms A/A)",
  `7` = "post A/B 25 ms",
  `8` = "baseline 5 (post 25 ms A/B)",
  `9` = "baseline 6 (post 25 ms A/B)",
  `10` = "post A/A 25 ms"
  )

blockType = c('baseline','baseline','A/A 200','baseline',
              'A/B 200','baseline','A/B 25','baseline','baseline','A/B 25')

whichCompareVec = list(c(2,3),c(4,5),c(6,7),c(9,10))
