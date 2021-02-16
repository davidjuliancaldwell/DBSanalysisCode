subjectNum = 19
chanIntVec = c(5,8)
chanIntConditioningPairVec = c(1,0)
blockIntLM = c(1,2,3,4,5,6,7,8)

blockIntPlot = blockIntLM

blockNames <- c(
  `1` = "baseline 1",
  `2` = "baseline 2 (pre conditioning)",
  `3` = "post A/B 100 ms",
  `4` = "baseline 3 (post 100 ms A/B)",
  `5` = "post A/A 100 ms",
  `6` = "baseline 4 (post 100 ms A/A)",
  `7` = "post A/B 200 ms",
  `8` = "baseline 5 (post 200 ms A/B)",
  `9` = "during DBS",
  `10`= "post DBS 1",
  `11` = "post DBS 2"
  )

#blockType = c('baseline','baseline','A/B 100','baseline',
#              'A/A 100','baseline','A/B 200','baseline',
 #             'during DBS','post DBS','baseline')


blockType = c('baseline','baseline','A/B 100','baseline',
              'A/A 100','baseline','A/B 200','baseline')

whichCompareVec = list(c(2,3),c(4,5),c(6,7))
