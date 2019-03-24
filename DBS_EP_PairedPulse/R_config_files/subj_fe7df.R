subjectNum = 7
chanIntVec = c(6)
chanIntConditioningPairVec = c(1)
blockIntLM = c(1,2,3,4,5,6,7,8,9)

blockIntPlot = blockIntLM
blockNames <- c(
  `1` = "baseline 1",
  `2` = "baseline 2",
  `3` = "baseline 3",
  `4` = "post 25 ms A/B 1",
  `5` = "baseline 4 post 25 ms A/B 2",
  `6` = "post 200 ms A/B 1",
  `7` = "baseline 5 post 200 ms A/B 2",
  `8` = "baseline 6 post 200 ms A/B 3",
  `9` = "baseline 7 post 200 ms A/B 4"
)

blockType = c('baseline','baseline','baseline','A/B 25','baseline',
              'A/B 200','baseline','baseline','baseline')

whichCompareVec = list(c(3,4),c(5,6))
