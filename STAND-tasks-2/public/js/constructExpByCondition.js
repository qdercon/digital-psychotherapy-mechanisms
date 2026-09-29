// construct experiment from various jsPsych and phaser components

// import task order randomization
import { taskCond } from "./versionInfo.js"; 

// import different task run functions
import { runStudyRew } from "./rew-eff/task.js";
import { runStudyCaus } from "./causal-attr/constructStudy.js";

// construct study depending on task order randomisation (all participants will do both tasks, we just
// want to randomise which they do first in case of order effects):
export function runStudy(){

	// to do the causal attribution task first 
	if ( taskCond == "causal-attr-first" ) {
		// start in jspsych div for causal attr instructions, task
		document.getElementById('game-container').style.display = "none";     // hide phaser container
		document.body.style.background = "aliceblue";
		runStudyCaus();
	}
	
	//
	if ( taskCond == "rew-eff-first" ) {
		// start in phaser div for instructions, practice, task 
		document.getElementById('game-container').style.display = "flex";     // show phaser container
		document.getElementById('jspsych-target').style.display = "none";     // hide jspsych container
		runStudyRew();
	}
}


