subjectNum = 8
chanIntVec = c(6)
chanIntConditioningPairVec = c(1)
blockIntLM = c(6,7,9,11,13,15)

blockIntPlot = blockIntLM

blockNames <- c(
  `6` = "baseline 1",
  `7` = "baseline 2 (pre conditioning)",
  `9` = "post A/B 200 ms",
  `11` = "baseline 3 (post 200 ms A/B)",
  `13` = "post A/A 25 ms",
  `15` = "post A/B 25 ms"
)

blockType = c('baseline','baseline','A/B 200','baseline',
              'A/A 25','A/B 25')

whichCompareVec = list(c(7,9),c(11,13),c(11,15))
