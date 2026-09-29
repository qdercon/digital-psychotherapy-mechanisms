// set up variables describing this specific task version

// task version
var version = "STAND-tasks-1";			       // experiment version (used to create data collection in firestore)
var infoSheet = "../assets/participant-information-sheet-230209-STAND.pdf";  
var briefStudyDescr = "For this part of the study (<b>part one</b>), we will ask you to complete <b>two different games</b>. "+
					  "At the end, you will also be also asked to <b>answer some questions about yourself, your feelings, and your mood</b>. ";

// are we debugging, or running for real?
var debugging = false;							// !!set to "false" for real exp!!
let allowDevices = false;                		// allow participants to access this task on mobile devices?

// time and payment info for this task version
var approxTime = 45;   			 				// approx time to complete this version of the experiment (minutes)
var hourlyRate = 7.5;							// 7.50 hourly rate (GBP)
var baseEarn = ((approxTime/60)*hourlyRate);    
var nQuests = 6;  								// how many questionnaires will we ask participant to complete?      

// set choice task variables
var nBlocksChoice = 2;

// set effort-related task variables
var nBlocksRew = 4;
var effortTime = 10000;					   		// time participant will have to try and exert effort (ms)
var pracTrialEfforts = [56, 26, 52, 68, 53];    // array of efforts ppts will be asked to perform in effortTime during practice
var minPressMax = 55;        				    // set a minimum on max press count to avoid gaming the practice trials
var bonusRate = 0.2;			 	 			// additional bonus per coin collected (GBPpence)
const maxCoins = 266;
var maxBonus = (maxCoins*bonusRate*2)/100;
var approxTimeRewTask = 12;                     // approx time to complete each rew-eff task

// randomly assign task order condition
var taskCond;
let r1 = Math.random();                          // (no random seed availble for this rng as this depends on browser)
if (r1 < 0.5) {
	taskCond = "rew-eff-first";
} else {
	taskCond = "causal-attr-first";
}

// randomly assign intervention (STAND module) condition
var intCond;
let r2 = Math.random();                          // (no random seed availble for this rng as this depends on browser)
if (r2 < 0.5) {
	intCond = "BA";
} else {
	intCond = "CR";
}

export { version, infoSheet, briefStudyDescr, debugging, allowDevices, 
		 approxTime, hourlyRate, baseEarn, bonusRate, nQuests, 
		 nBlocksChoice, nBlocksRew, effortTime, pracTrialEfforts, minPressMax, maxBonus, approxTimeRewTask,
		 taskCond, intCond
		};