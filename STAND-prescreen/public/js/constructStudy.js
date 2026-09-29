// import task info from versionInfo file
import { allowDevices } from "./versionInfo.js"; 
 
// import study elements
import { timeline_quests } from "./selfReports.js";

////////////////////////////////////// initialise jsPsych ///////////////////////////////////
var jsPsych = initJsPsych({
    display_element: 'jspsych-target',
    // show_progress_bar: false  
    show_progress_bar: true,      
    message_progress_bar: 'progress',
    auto_update_progress_bar: false
});
// export jsPsych object so can be accessed by other study modules
export { jsPsych };

//////////////////////////// construct overall study timeline ///////////////////////////////
export function runStudy(){
    // initialise overall study timeline
    var timeline = [];

    // construct the study timeline...           
    // 1. instructions and self-reports, debrief and study end screens
    timeline = timeline.concat(timeline_quests);                    
    
    // ...and run it!
    jsPsych.run(timeline);

}; 
