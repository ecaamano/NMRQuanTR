
# Tests for functions in lib.R

data = read.csv('test_data/dataForQuant/test_data_1', header=T)
specNom = read.csv('test_data/SpectraNomenclature.csv', header=T, stringsAsFactors = F)
specNom_CC = read.csv('test_data/TemplateCalCurves.csv', header=T, stringsAsFactors = F)
carbToMet = read.csv('test_data/carbonToMetabolite.csv', header=T, stringsAsFactors = F)
calcurveCoeffs = read.csv('test_data/calCurveCoeffs.csv', header=T, stringsAsFactors = F)

data_1 <- parseSpinSysData(data)
data_2 <- changeSpins2Names(data_1, specNom)
data_3 <- normaliseByPeak(data_2, normBy = 'TSP', exclude = 'none', concRef = 2.5, concRefSC = 2.5)
data_4 <- changeNomenclature(data_3, specNom, excludeNotFound=T)
data_5 <- setNamesAndOrdering(data_4)
data_6 <- calcConc(data_5, calcurveCoeffs)
data_7 <- adjustAliquot(data_6, vNMRtube = 300, vAliquot = 800)
  
test_that("Test Data", {
  expect_equal(ncol(data),8)
  expect_equal(nrow(data),39)
})

test_that("Test SpecNom", {
  expect_equal(ncol(specNom),5)
  expect_equal(nrow(specNom),95)
  expect_type(specNom[,1], 'character')
  expect_type(specNom[,2], 'integer')
})

test_that("Test parseSpinSysData", {
  expect_equal(ncol(data_1),9)
  expect_equal(nrow(data_1),nrow(data))
  expect_equal(names(data_1)[5],'Sp.system')
  expect_equal(names(data_1)[6],'Assign.F1')
  expect_equal(names(data_1)[7],'Assign.F2')
  expect_type(data_1, 'list')
})

test_that("Test changeSpins2Names", {
  expect_type(data_2$Sp.system, 'character')
  expect_equal(ncol(data_1), ncol(data_2))
  expect_true(all(names(data_2)==names(data_2)))
})

test_that("Test NormaliseByPeak", {
  expect_equal(nrow(data_3), nrow(data_2)-1)
  expect_equal(ncol(data_3), ncol(data_2))
  expect_true(all(names(data_3)==names(data_2)))
})

test_that("Test changeNomenclature", {
  expect_equal(ncol(data_4), 3)
  expect_true(all(names(data_4)==c('Carbon', 'Height', 'Volume')))
})

test_that("Test setNamesAndOrdering", {
  expect_equal(nrow(data_5), nrow(data_4))
  expect_equal(ncol(data_5), ncol(data_4))
  expect_true(all(names(data_5)==names(data_4)))
})

test_that("Test calcConc", {
  expect_equal(nrow(data_6), nrow(data_5))
  expect_equal(ncol(data_6), ncol(data_5))
  expect_true(all(names(data_6)==c('Carbon', 'Conc_Height', 'Conc_Volume')))
})

test_that("Test adjustAliquot", {
  expect_equal(nrow(data_7), nrow(data_6))
  expect_equal(ncol(data_7), ncol(data_6))
  expect_true(all(names(data_7)==names(data_6)))
})

test_that("Test ", {
  expect_equal(nrow(data_1), nrow(data_2))
  expect_equal(ncol(data_1), ncol(data_2))
  expect_true(all(names(data_2)==names(data_2)))
})