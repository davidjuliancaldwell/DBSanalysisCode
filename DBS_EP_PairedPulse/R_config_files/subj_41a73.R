subjectNum = 14
chanIntVec = c(5)
chanIntConditioningPairVec = c(1)
blockIntLM = c(1,2,3,4,5,6,7,8,9,10)

blockIntPlot = blockIntLM

blockNames <- c(
  `1` = "baseline 1",
  `2` = "baseline 2 (pre conditioning)",
  `3` = "post A/B 200 ms",
  `4` = "baseline 3 (post A/B 200)",
  `5` = "baseline 4",
  `6` = "post A/A 200 ms",
  `7` = "baseline 5 (post A/A 200)",
  `8` = "baseline 6",
  `9` = "baseline 7",
  `10` = "baseline 8"
  )

blockType = c('baseline','baseline','A/B 200','baseline',
              'baseline','A/A 200','baseline','baseline','baseline','baseline')

whichCompareVec = list(c(2,3),c(5,6))
