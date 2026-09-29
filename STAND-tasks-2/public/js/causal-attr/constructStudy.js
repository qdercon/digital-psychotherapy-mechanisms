// import study elements
import { timeline_instructions_choice } from "./instructionsChoice.js";
import { timeline_choice_2 } from "./taskChoice.js";

////////////////////////////////////// initialise jsPsych ///////////////////////////////////
var jsPsych = initJsPsych({
    display_element: 'jspsych-target',
    show_progress_bar: false  
    // // show_progress_bar: true,      
    // message_progress_bar: 'progress',
    // auto_update_progress_bar: false
});
// export jsPsych object so can be accessed by other study modules
export { jsPsych };

//////////////////////////// construct overall CA first timeline ///////////////////////////////
export function runStudyCaus(){
    // initialise overall study timeline
    var timeline = [];

    // construct the study timeline...
    // initial instructions, choice test 1
    timeline = timeline.concat(timeline_instructions_choice);   //  timeline_instructions_choice_1
    timeline = timeline.concat(timeline_choice_2);                                              
    
    // ...and run it!
    jsPsych.run(timeline);

}; 