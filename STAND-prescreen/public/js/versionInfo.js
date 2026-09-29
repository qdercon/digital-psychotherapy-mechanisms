// set up variables describing this specific task version

// task version
var version = "STAND-prescreen";			    // experiment version (used to create data collection in firestore)
var infoSheet = "../assets/participant-information-sheet-230209-noSTAND.pdf";  
var briefStudyDescr = "For this particular study, we will ask you to answer some questions about your feelings and mood.";

// are we debugging, or running for real?
var debugging = false;							// !!set to "false" for real exp!!
let allowDevices = false;                		// allow participants to access this task on mobile devices?

// time and payment info for this task version
var approxTime = 5;   			 				// approx time to complete this version of the experiment (minutes)
var hourlyRate = 9.0;							// 7.50 hourly rate (GBP)
var baseEarn = ((approxTime/60)*hourlyRate);    
var nQuests = 3;  								// how many questionnaires will we ask participant to complete?      


export { version, infoSheet, briefStudyDescr, debugging, allowDevices, 
		 approxTime, hourlyRate, baseEarn, nQuests
		};