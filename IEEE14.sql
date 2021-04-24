/*
 Navicat Premium Data Transfer

 Source Server         : Server
 Source Server Type    : MySQL
 Source Server Version : 80022
 Source Host           : wadegao.tpddns.net:3306
 Source Schema         : IEEE14

 Target Server Type    : MySQL
 Target Server Version : 80022
 File Encoding         : 65001

 Date: 24/04/2021 21:47:21
*/

SET NAMES utf8mb4;
SET FOREIGN_KEY_CHECKS = 0;

-- ----------------------------
-- Table structure for branch
-- ----------------------------
DROP TABLE IF EXISTS `branch`;
CREATE TABLE `branch` (
  `tapBusNo` int unsigned DEFAULT NULL,
  `ZbusNo` int unsigned DEFAULT NULL,
  `loadFlowAreaNo` int unsigned DEFAULT NULL,
  `lossZoneNo` int unsigned DEFAULT NULL,
  `circuit` int unsigned DEFAULT NULL,
  `type` tinyint unsigned DEFAULT NULL,
  `R` float(15,10) DEFAULT NULL,
  `X` float(15,10) DEFAULT NULL,
  `B` float(15,10) DEFAULT NULL,
  `lineMVAratingNo1` float(15,10) DEFAULT NULL,
  `lineMVAratingNo2` float(15,10) DEFAULT NULL,
  `lineMVAratingNo3` float(15,10) DEFAULT NULL,
  `ctrlBusNo` int unsigned DEFAULT NULL,
  `side` tinyint unsigned DEFAULT NULL,
  `transformerFinalTurnsRatio` float(15,10) DEFAULT NULL,
  `transformerPhaseShifterFinalAngle` float(15,10) DEFAULT NULL,
  `minTapOrPhaseShift` float(15,10) DEFAULT NULL,
  `maxTapOrPhaseShift` float(15,10) DEFAULT NULL,
  `stepSize` float(15,10) DEFAULT NULL,
  `minVolt` float(15,10) DEFAULT NULL,
  `maxVolt` float(15,10) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

-- ----------------------------
-- Records of branch
-- ----------------------------
BEGIN;
INSERT INTO `branch` VALUES (1, 2, 1, 1, 1, 0, 0.0193799995, 0.0591700003, 0.0527999997, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (1, 5, 1, 1, 1, 0, 0.0540300012, 0.2230399996, 0.0491999984, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (2, 3, 1, 1, 1, 0, 0.0469899997, 0.1979700029, 0.0438000001, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (2, 4, 1, 1, 1, 0, 0.0581099987, 0.1763200015, 0.0340000018, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (2, 5, 1, 1, 1, 0, 0.0569499992, 0.1738799959, 0.0346000008, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (3, 4, 1, 1, 1, 0, 0.0670100003, 0.1710299999, 0.0127999997, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (4, 5, 1, 1, 1, 0, 0.0133499997, 0.0421099998, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (4, 7, 1, 1, 1, 1, 0.0000000000, 0.2091200054, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.9779999852, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (4, 9, 1, 1, 1, 1, 0.0000000000, 0.5561800003, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.9689999819, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (5, 6, 1, 1, 1, 1, 0.0000000000, 0.2520200014, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.9319999814, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (6, 11, 1, 1, 1, 0, 0.0949800014, 0.1988999993, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (6, 12, 1, 1, 1, 0, 0.1229100004, 0.2558099926, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (6, 13, 1, 1, 1, 0, 0.0661500022, 0.1302700043, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (7, 8, 1, 1, 1, 1, 0.0000000000, 0.1761499941, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (7, 9, 1, 1, 1, 1, 0.0000000000, 0.1100099981, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (9, 10, 1, 1, 1, 0, 0.0318100005, 0.0844999999, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (9, 14, 1, 1, 1, 0, 0.1271100044, 0.2703799903, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (10, 11, 1, 1, 1, 0, 0.0820500031, 0.1920700073, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (12, 13, 1, 1, 1, 0, 0.2209199965, 0.1998800039, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
INSERT INTO `branch` VALUES (13, 14, 1, 1, 1, 0, 0.1709299982, 0.3480199873, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0, 0, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000);
COMMIT;

-- ----------------------------
-- Table structure for bus
-- ----------------------------
DROP TABLE IF EXISTS `bus`;
CREATE TABLE `bus` (
  `bus` int unsigned DEFAULT NULL,
  `name` varchar(20) DEFAULT NULL,
  `loadFlowAreaNo` int unsigned DEFAULT NULL,
  `lossZoneNo` int unsigned DEFAULT NULL,
  `type` tinyint unsigned DEFAULT NULL,
  `finalVolt` float(15,10) DEFAULT NULL,
  `finalAngle` float(15,10) DEFAULT NULL,
  `loadMW` float(15,10) DEFAULT NULL,
  `loadMVAR` float(15,10) DEFAULT NULL,
  `generationMW` float(15,10) DEFAULT NULL,
  `generationMVAR` float(15,10) DEFAULT NULL,
  `baseKV` float(15,10) DEFAULT NULL,
  `desiredVolt` float(15,10) DEFAULT NULL,
  `maxMVAR` float(15,10) DEFAULT NULL,
  `minMVAR` float(15,10) DEFAULT NULL,
  `G` float(15,10) DEFAULT NULL,
  `B` float(15,10) DEFAULT NULL,
  `remoteCtrlBusNo` int unsigned DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

-- ----------------------------
-- Records of bus
-- ----------------------------
BEGIN;
INSERT INTO `bus` VALUES (1, 'Bus 1', 1, 1, 3, 1.0599999428, 0.0000000000, 0.0000000000, 0.0000000000, 232.3999938965, -16.8999996185, 132.0000000000, 1.0599999428, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (2, 'Bus 2', 1, 1, 2, 1.0449999571, -4.9800000191, 21.7000007629, 12.6999998093, 40.0000000000, 42.4000015259, 132.0000000000, 1.0449999571, 50.0000000000, -40.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (3, 'Bus 3', 1, 1, 2, 1.0099999905, -12.7200002670, 94.1999969482, 19.0000000000, 0.0000000000, 23.3999996185, 132.0000000000, 1.0099999905, 40.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (4, 'Bus 4', 1, 1, 0, 1.0190000534, -10.3299999237, 47.7999992371, -3.9000000954, 0.0000000000, 0.0000000000, 132.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (5, 'Bus 5', 1, 1, 0, 1.0199999809, -8.7799997330, 7.5999999046, 1.6000000238, 0.0000000000, 0.0000000000, 132.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (6, 'Bus 6', 1, 1, 2, 1.0700000525, -14.2200002670, 11.1999998093, 7.5000000000, 0.0000000000, 12.1999998093, 35.0000000000, 1.0700000525, 24.0000000000, -6.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (7, 'Bus 7', 1, 1, 0, 1.0620000362, -13.3699998856, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 35.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (8, 'Bus 8', 1, 1, 2, 1.0900000334, -13.3599996567, 0.0000000000, 0.0000000000, 0.0000000000, 17.3999996185, 10.0000000000, 1.0900000334, 24.0000000000, -6.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (9, 'Bus 9', 1, 1, 0, 1.0559999943, -14.9399995804, 29.5000000000, 16.6000003815, 0.0000000000, 0.0000000000, 35.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.1899999976, 0);
INSERT INTO `bus` VALUES (10, 'Bus 10', 1, 1, 0, 1.0509999990, -15.1000003815, 9.0000000000, 5.8000001907, 0.0000000000, 0.0000000000, 35.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (11, 'Bus 11', 1, 1, 0, 1.0570000410, -14.7899999619, 3.5000000000, 1.7999999523, 0.0000000000, 0.0000000000, 35.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (12, 'Bus 12', 1, 1, 0, 1.0549999475, -15.0699996948, 6.0999999046, 1.6000000238, 0.0000000000, 0.0000000000, 35.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (13, 'Bus 13', 1, 1, 0, 1.0499999523, -15.1599998474, 13.5000000000, 5.8000001907, 0.0000000000, 0.0000000000, 35.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0);
INSERT INTO `bus` VALUES (14, 'Bus 14', 1, 1, 0, 1.0360000134, -16.0400009155, 14.8999996185, 5.0000000000, 0.0000000000, 0.0000000000, 35.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0);
COMMIT;

SET FOREIGN_KEY_CHECKS = 1;
