subjectNum = 9
chanIntVec = c(4)
blockIntLM = c(2,3,5,7,11)
blockIntLM = c(2,3,4,5,6,7,10,11)

blockIntPlot = c(2,3,4,5,6,7,10,11)
blockNames <- c(
  `2` = "baseline 2 (pre conditioning)",
  `3` = "post A/B 25 ms",
  `4` = "baseline 3 (post 25 ms)",
  `5` = "post A/A 25 ms",
  `6` = "baseline 4 (post 25 ms A/A)",
  `7` = "post A/B 200 ms",
  `10` = "'post A/B 200 ms 12 minutes later",
  `11` = "post A/B 25 ms second time"
  )

blockType = c('baseline','A/B 25','baseline',
              'A/A 25','baseline','A/B 200','A/B 200','A/B 25')

whichCompareVec = list(c(2,3),c(4,5),c(6,7),c(6,10,11))
