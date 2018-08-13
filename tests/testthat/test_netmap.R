context("netmap")

test_that("npn map construction method works ", 
  {
	 data(CviCol)
	 M1 <- netmap(CviCol, method="npn", cross= "inbred", rho = 0.4)
	 #dput(M1$map)
	
	  expect_equal( M1$map, structure(list(markers = structure(c(1L, 2L, 3L, 4L, 5L, 6L, 
		7L, 8L, 9L, 11L, 10L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 
		20L, 21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 30L, 31L, 32L, 
		33L, 34L, 35L, 36L, 37L, 38L, 39L, 40L, 41L, 42L, 43L, 44L, 45L, 
		46L, 47L, 48L, 49L, 50L, 51L, 52L, 53L, 54L, 55L, 56L, 57L, 58L, 
		59L, 60L, 61L, 62L, 63L, 64L, 65L, 66L, 67L, 68L, 69L, 70L, 71L, 
		72L, 73L, 74L, 75L, 76L, 77L, 78L, 79L, 80L, 81L, 82L, 83L, 84L, 
		85L, 86L, 87L, 88L, 89L, 90L), .Label = c("c1_00593", "c1_02212", 
		"c1_02992", "c1_04176", "c1_05593", "c1_08385", "c1_09782", "c1_11160", 
		"c1_12295", "c1_13869", "c1_13926", "c1_15634", "c1_16875", "c1_18433", 
		"c1_19478", "c1_20384", "c1_22181", "c1_23381", "c1_24795", "c1_25698", 
		"c1_26993", "c1_28454", "c1_28667", "c1_29898", "c2_00593", "c2_02365", 
		"c2_03041", "c2_04263", "c2_06280", "c2_07650", "c2_10250", "c2_11457", 
		"c2_12435", "c2_13472", "c2_15252", "c2_16837", "c2_17606", "c2_18753", 
		"c3_00580", "c3_00885", "c3_01901", "c3_02968", "c3_04141", "c3_05141", 
		"c3_06631", "c3_08042", "c3_09748", "c3_10996", "c3_11192", "c3_12647", 
		"c3_15117", "c3_16677", "c3_18180", "c3_20729", "c3_22147", "c4_00012", 
		"c4_00641", "c4_02133", "c4_03833", "c4_04877", "c4_05850", "c4_06923", 
		"c4_07549", "c4_07740", "c4_08930", "c4_10609", "c4_11878", "c4_13171", 
		"c4_14819", "c4_15765", "c4_17684", "c5_00576", "c5_01587", "c5_02900", 
		"c5_04011", "c5_05319", "c5_06820", "c5_07442", "c5_08563", "c5_10428", 
		"c5_13614", "c5_14766", "c5_17570", "c5_19316", "c5_20318", "c5_21319", 
		"c5_22415", "c5_23116", "c5_24997", "c5_26671"), class = "factor"), 
			LG = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
			1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 
			2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 
			3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 
			4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 
			5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 
			5L, 5L)), .Names = c("markers", "LG"), row.names = c(NA, 
		-90L), class = "data.frame")  )  
   })


test_that("Gibbs map construction method works ", 
  {
	 data(CviCol)
	 D2 <- CviCol[ 1:20, c(1:5, 61:65)]
	 M2 <- netmap(D2, method="gibbs", cross= "inbred", rho = 0.3, ncores = 1)
	
	 expect_equal( max(M2$map$LG), 2 )  
   })


test_that("Approx map construction method works ", 
  {
	 data(CviCol)
	 D3 <- CviCol[ 1:5, 1:5 ]
	 M3 <- netmap(D3, method="approx", cross= "inbred", rho = 0.6, ncores = 1)
	
	 expect_equal( max(M3$map$LG), 2 )  
   })








