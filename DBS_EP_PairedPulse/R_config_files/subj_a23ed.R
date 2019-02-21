subjectNum = 20
chanIntVec = c(5,8)
chanIntConditioningPairVec = c(1,0)
blockIntLM = c(1,2,3,4,5,6,7,8,9,10,11)

blockIntPlot = blockIntLM

blockNames <- c(
  `1` = "baseline 1",
  `2` = "baseline 2 (pre conditioning)",
  `3` = "post A/A 100 ms",
  `4` = "baseline 3 (post 100 ms A/A)",
  `5` = "post A/B 100 ms",
  `6` = "baseline 4 (post 100 ms A/B)",
  `7` = "post A/A 200 ms",
  `8` = "baseline 5 (post 200 ms A/A)",
  `9` = "post A/B 200 ms",
  `10`= "baseline 6 (post 200 ms A/B)",
  `11` = "baseline 7 - pre DBS"
  )

blockType = c('baseline','baseline','A/A 100','baseline',
              'A/B 100','baseline','A/A 200','baseline','A/B 200','baseline','baseline')

whichCompareVec = list(c(2,3),c(4,5),c(6,7),c(9,10))
