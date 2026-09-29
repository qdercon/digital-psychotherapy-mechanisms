// Script to run a series of self-report questionnaires, using built-in JsPsych functionality

// import relevant task info
import { nQuests} from "./versionInfo.js";

// import data saving functions
import { saveQuestData, saveEndData } from "./saveData.js";

// import global jspsych object
import { jsPsych } from "./causal-attr/constructStudy.js";

// initialize sizing vars
var scaleDisplayWidth = 600;  // in px

///////////////////////////////////////////// MISC TEXT /////////////////////////////////////////////////////////
var questsIntroText = {
  type: jsPsychHtmlButtonResponse,
  choices: ['start'],
  is_html: true,
  stimulus: ("<p><h2>Thank you!</h2></p>"+
            "<br>"+
            "<p>"+
            "<b>For the final part of the study, we woud like you to answer some questions about "+
            "yourself, your feelings and mood</b>. "+
            "</p>"+
            "<p>"+
            "We appreciate that you have already done a lot of activities at this point, "+
            "but it is <i>really important for the aims of our study that you answer each of the remaining questions as "+
            "accurately and truthfully as possible</i>. This will allow us to learn as much as possible "+
            "from the data you have already provided during the other parts of the study. "+
            "</p>"+
            "<p>"+
            "In total, we will ask you to complete <b>"+nQuests+" short questionnaires (all 10 "+
            "questions or less), and some questions about yourself and your personal circumstances</b>."+
            "</p>"+
            "<p>"+
            // "<b>The 'progress' bar at the top of the page shows you how many of the questions you have "+
            // "left to go.</b>"+
            "<br><br><br><br><br><br>"+
            "</p>"),
  on_start: function() {
    //this.type.jsPsych.setProgressBar(0);
    document.body.style.background = "aliceblue";
  }
};

var questsEndScreen = {
  type: jsPsychHtmlButtonResponse,
  timing_post_trial: 0,
  choices: ['finish the study'],
  is_html: true,
  stimulus: ("<p>"+
            "<h2>Thank you very much for your time. You have now finished the study.</h2>"+
            // "</p>"+
            // "<b>You will receive instructions via Prolific letting you know how to access the next "+
            // "part of the study</b> (the online course)."+
            // "</p>"+
            "</p>"+
            "If you would like to find out more about the ideas behind this study, "+
            "please see <a href=\"https://www.psychologytools.com/self-help/thoughts-in-cbt/\" target=\"_blank\">this article</a> "+
            "about why some psychologists believe the way we <i>interpret</i> events is key to understanding our feelings about them, "+
            "and <a href=\"https://www.psychologytools.com/self-help/behavioral-activation/\" target=\"_blank\">this article</a> about why some psychologist "+
            "think that changing the way we <i>act</i> can help improve our mood. "+
            "</p>"+
            "<p>"+
            "If you became upset at any point during the study, "+
            "or are concerned about your mental health for any other reason, we recommend the below resources for further "+
            "information. You may also wish to discuss any concerns with your family doctor."+ 
            "</p>  " +
            "<ul>  " +
            "<p><li><a href=\"http://mind.org.uk\" target=\"_blank\">Mind Charity</a></li></p>"+
            "<p><li><a href=\"https://www.samaritans.org\" target=\"_blank\">The Samaritans</a></li></p>"+
            "<p><li><a href=\"https://www.nhs.uk/mental-health\" target=\"_blank\">NHS Choices mental health page</a></li></p>"+
            "</ul>"+
            "</p>"+
            "<b>Please click the button below to submit your data back to Prolific!</b>"+
            // ! use target=blank to ensure links open in new window, and don't mess with prolific submission
            "</p>"),
  on_finish: function() { 
    // redirect to study completion page 
    window.location = "https://app.prolific.com/submissions/complete?cc=C1MC610R";
  } 
};  

///////////////////////////////////////////// SCALES /////////////////////////////////////////////////////////
//////////////////////////// DAS-SF ///////////////////////
// define response labels
var respOptsDAS = ["Totally Agree", "Agree", "Disagree", "Totally Disagree"];

// scale items
var DAS = {
  type: jsPsychSurveyLikert,
  preamble: ("<p>"+
            "The sentences below describe people’s attitudes. Please select how much each sentence describes "+
             "your attitude. Your answer should describe the way you think <b>most of the time</b>."+
             "</p>"),
  questions: [
    {prompt: "<b>If I don’t set the highest standards for myself, I am likely to end up a second-rate person</b>", 
      name: "DAS_1", labels: respOptsDAS, required:true, horizontal: true}, 
    {prompt: "<b>My value as a person depends greatly on what others think of me</b>", 
      name: "DAS_2", labels: respOptsDAS, required: true, horizontal: true},
    {prompt: "<b>People will probably think less of me if I make a mistake</b>", 
      name: "DAS_3", labels: respOptsDAS, required:true, horizontal: true},
    {prompt: "<b>I am nothing if a person I love doesn’t love me</b>", 
      name: "DAS_4", labels: respOptsDAS, required:true, horizontal: true},
    {prompt: "<b>I think it should be against the law to listen to music.</b>", 
      name: "catch_2", labels: respOptsDAS, required:true, horizontal: true},
    {prompt: "<b>If other people know what you are really like, they will think less of you</b>", 
      name: "DAS_5", labels: respOptsDAS, required:true, horizontal: true},
    {prompt: "<b>If I fail at my work, then I am a failure as a person</b>", 
      name: "DAS_6", labels: respOptsDAS, required:true, horizontal: true},
    {prompt: "<b>My happiness depends more on other people than it does me</b>", 
      name: "DAS_7", labels: respOptsDAS, required:true, horizontal: true},
    {prompt: "<b>I cannot be happy unless most people I know admire me.</b>", 
      name: "DAS_8", labels: respOptsDAS, required:true, horizontal: true},
    {prompt: "<b>It is best to give up your own interests in order to please other people.</b>", 
      name: "DAS_9", labels: respOptsDAS, required:true, horizontal: true}
  ],
  button_label: 'continue',
  scale_width: scaleDisplayWidth,  
  on_finish: function () {
    // set progress bar manually
    // this.type.jsPsych.setProgressBar(0.4);
    // get response and RT data
    var respData = this.type.jsPsych.data.getLastTrialData().trials[0].response;
    var respRT = this.type.jsPsych.data.getLastTrialData().trials[0].rt;
    saveQuestData("DAS", respData, respRT);
  }
};

//////////////////////////// PHQs ///////////////////////
// define response labels
var respOptsPHQ9 = ["Not at all", "Several days", "More than half the days", "Nearly every day"];
var respOptsPHQ9Difficulty = ["Not difficult at all", "Somewhat difficult", "Very difficult", "Extremely difficult"];

// PHQ9 scale items
var PHQ9 = {
  type: jsPsychSurveyLikert,
  preamble: ("Over the <b>last week</b>, how often have you been bothered by any of the following problems?\n "),
  questions: [
    {prompt: "<b>Little interest or pleasure in doing things</b>", 
      name: "PHQ9_1", labels: respOptsPHQ9, required:true, horizontal: true}, 
    {prompt: "<b>Feeling down, depressed, or hopeless</b>", 
      name: "PHQ9_2", labels: respOptsPHQ9, required: true, horizontal: true},
    {prompt: "<b>Trouble falling/staying asleep, sleeping too much</b>", 
      name: "PHQ9_3", labels: respOptsPHQ9, required:true, horizontal: true},
    {prompt: "<b>Feeling tired or having little energy</b>", 
      name: "PHQ9_4", labels: respOptsPHQ9, required:true, horizontal: true},
    {prompt: "<b>Poor appetite or overeating</b>", 
      name: "PHQ9_5", labels: respOptsPHQ9, required:true, horizontal: true},
    {prompt: "<b>Feeling bad about yourself or that you are a failure or have let yourself or your family down</b>", 
      name: "PHQ9_6", labels: respOptsPHQ9, required:true, horizontal: true},
    {prompt: "<b>Trouble concentrating on things, such as reading the newspaper or watching television.</b>", 
      name: "PHQ9_7", labels: respOptsPHQ9, required:true, horizontal: true},
    {prompt: "<b>Moving or speaking so slowly that other people could have noticed.\n"+
            "Or the opposite; being so fidgety or restless that you have been moving around a lot more than usual.</b>", 
      name: "PHQ9_8", labels: respOptsPHQ9, required:true, horizontal: true},
    {prompt: "<b>Thoughts that you would be better off dead or of hurting yourself in some way.</b>", 
      name: "PHQ9_9", labels: respOptsPHQ9, required:true, horizontal: true},
    {prompt: "If you have been bothered by any of the above, how difficult have these problems "+
             "made it for you to do your work, take care of things at home, or get along with other people?", 
      name: "PHQ9_D", labels: respOptsPHQ9Difficulty, required:true, horizontal: true}
  ],
  button_label: 'continue',
  scale_width: scaleDisplayWidth,  
  on_finish: function(){
    // get response and RT data 
    var respData = this.type.jsPsych.data.getLastTrialData().trials[0].response;
    var respRT = this.type.jsPsych.data.getLastTrialData().trials[0].rt;
    saveQuestData("PHQ9", respData, respRT);
  }
};

////////////////////MINI-SPIN//////////////////////////////
// define response labels
var respOptsMiniSPIN = ["Not at all", "A little bit", "Somewhat", "Very much", "Extremely"];

// scale items
var miniSPIN = {
  type: jsPsychSurveyLikert,
  preamble: ("<p>"+
             "For each statement below, please select how well it describes you. "+
             "</p>"),
  questions: [
    {prompt: "<b>Fear of embarrassment causes me to avoid doing things or speaking to people.</b>", 
    name: "miniSPIN_1", labels: respOptsMiniSPIN, required:true, horizontal: true}, 
    {prompt: "<b>I avoid activities in which I am the center of attention.</b>", 
    name: "miniSPIN_2", labels: respOptsMiniSPIN, required:true, horizontal: true}, 
    {prompt: "<b>Being embarrassed or looking stupid are among my worst fears.</b>", 
    name: "miniSPIN_3", labels: respOptsMiniSPIN, required:true, horizontal: true}
  ],
  button_label: 'continue',
  scale_width: scaleDisplayWidth,
  on_finish: function() {
    // // set progress bar manually
    // this.type.jsPsych.setProgressBar(0.8);
    // get response and RT data
    var respData = this.type.jsPsych.data.getLastTrialData().trials[0].response;
    var respRT = this.type.jsPsych.data.getLastTrialData().trials[0].rt;
    saveQuestData("miniSPIN", respData, respRT);
  }
};

// //////////////////////////// SHAPS ///////////////////////
// // define response labels
// var respOptsSHAPS = ["Strongly disagree", "Disagree", "Agree", "Strongly agree"];

// // scale items
// var SHAPS = {
//   type: jsPsychSurveyLikert,
//   preamble: ("For each statement below, please select how much you agree or disagree, "+
//              "depending on <b>how you have generally felt over the past two weeks.</b>"),
//   questions: [
//     {prompt: "<b>I would enjoy my favourite television or radio programme.</b>", 
//     name: "SHAPS_1", labels: respOptsSHAPS, required:true, horizontal: true}, 
//     {prompt: "<b>I would enjoy being with my family or close friends.</b>", 
//     name: "SHAPS_2", labels: respOptsSHAPS, required: true, horizontal: true},
//     {prompt: "<b>I would find pleasure in my hobbies and pastimes.</b>", 
//     name: "SHAPS_3", labels: respOptsSHAPS, required:true, horizontal: true},
//     {prompt: "<b>I would be able to enjoy my favourite meal.</b>", 
//     name: "SHAPS_4", labels: respOptsSHAPS, required:true, horizontal: true},
//     {prompt: "<b>I would enjoy a warm bath or refreshing shower.</b>", 
//     name: "SHAPS_5", labels: respOptsSHAPS, required:true, horizontal: true,},
//     {prompt: "<b>I would find pleasure in the scent of flowers or the smell of a fresh sea breeze or freshly baked bread.</b>", 
//     name: "SHAPS_6", labels: respOptsSHAPS, required:true, horizontal: true},
//     {prompt: "<b>I would enjoy seeing other people's smiling faces.</b>", 
//     name: "SHAPS_7", labels: respOptsSHAPS, required:true, horizontal: true},
//     {prompt: "<b>I would enjoy looking smart when I have made an effort with my appearance.</b>", 
//     name: "SHAPS_8", labels: respOptsSHAPS, required:true, horizontal: true},
//     {prompt: "<b>I would enjoy reading a book, magazine or newspaper.</b>", 
//     name: "SHAPS_9", labels: respOptsSHAPS, required:true, horizontal: true},
//     // {prompt: "<b>I would consider myself to have taken part in the 1917 Summer Olympic Games.</b>", 
//     // name: "catch_2", labels: respOptsSHAPS, required:true, horizontal: true},
//     {prompt: "<b>I would enjoy a cup of tea or coffee or my favorite drink.</b>", 
//     name: "SHAPS_10", labels: respOptsSHAPS, required:true, horizontal: true},
//     {prompt: "<b>I would find pleasure in small things, e.g. bright sunny day, a telephone call from a friend.</b>", 
//     name: "SHAPS_11", labels: respOptsSHAPS, required:true, horizontal: true},
//     {prompt: "<b>I would be able to enjoy a beautiful landscape or view.</b>", 
//     name: "SHAPS_12", labels: respOptsSHAPS, required:true, horizontal: true},
//     {prompt: "<b>I would get pleasure from helping others.</b>", 
//     name: "SHAPS_13", labels: respOptsSHAPS, required:true, horizontal: true},
//     {prompt: "<b>I would feel pleasure when I receive praise from other people.</b>", 
//     name: "SHAPS_14", labels: respOptsSHAPS, required:true, horizontal: true}
//   ],
//   button_label: 'continue',
//   scale_width: scaleDisplayWidth,
//   on_finish: function(){ 
//     // get response and RT data
//     var respData = this.type.jsPsych.data.getLastTrialData().trials[0].response;
//     var respRT = this.type.jsPsych.data.getLastTrialData().trials[0].rt;
//     saveQuestData("SHAPS", respData, respRT);
//   }
// };

///////////////////////////////////////////// AMI /////////////////////////////////////////////////////////
// define response labels
var respOptsAMI = ["Completely UNTRUE", "Mostly untrue", "Neither true nor untrue", "Quite true", "Completely TRUE"];

// scale items
var AMI = {
  type: jsPsychSurveyLikert,
  preamble: ("<p>"+
             "For each statement below, please select how appropriately it describes you. "+
             "</p>"+
             "<p>"+
             "Select <i>Completely True</i> if the statement describes you perfectly, and <i>Completely "+
             "Untrue</i> if the statement does not describe you at all, "+
             "thinking about <b>the last week</b>. "+
             "</p>"),
  questions: [
    // {prompt: "<b>I feel sad or upset when I hear bad news.</b>", 
    // name: "AMI_1", labels: respOptsAMI, required:true, horizontal: true}, 
    // {prompt: "<b>I start conversations with random people.</b>", 
    // name: "AMI_2", labels: respOptsAMI, required: true, horizontal: true},
    // {prompt: "<b>I enjoy doing things with people I have just met.</b>", 
    // name: "AMI_3", labels: respOptsAMI, required:true, horizontal: true},
    // {prompt: "<b>I suggest activities for me and my friends to do.</b>", 
    // name: "AMI_4", labels: respOptsAMI, required:true, horizontal: true},
    {prompt: "<b>I make decisions firmly and without hesitation.</b>", 
    name: "AMI_5", labels: respOptsAMI, required:true, horizontal: true},
    // {prompt: "<b>After making a decision, I will wonder if I have made the wrong choice.</b>",
    // name: "AMI_6", labels: respOptsAMI, required:true, horizontal: true}, 
    // {prompt: "<b>I competed in the 1917 Summer Olympic Games.</b>", 
    // name: "catch_2", labels: respOptsAMI, required:true, horizontal: true},
    // {prompt: "<b>Based on the last two weeks, I would say I care deeply about how my loved ones think of me.</b>", 
    // name: "AMI_7", labels: respOptsAMI, required:true, horizontal: true},
    // {prompt: "<b>I go out with friends on a weekly basis.</b>", 
    // name: "AMI_8", labels: respOptsAMI, required:true, horizontal: true},
    {prompt: "<b>When I decide to do something, I am able to make an effort easily.</b>", 
    name: "AMI_9", labels: respOptsAMI, required:true, horizontal: true},
    {prompt: "<b>I don't like to laze around.</b>", 
    name: "AMI_10", labels: respOptsAMI, required:true, horizontal: true},
    {prompt: "<b>I get things done when they need to be done, without requiring reminders from others.</b>", 
    name: "AMI_11", labels: respOptsAMI, required:true, horizontal: true},
    {prompt: "<b>When I decide to do something, I am motivated to see it through to the end.</b>", 
    name: "AMI_12", labels: respOptsAMI, required:true, horizontal: true},
    // {prompt: "<b>I feel awful if I say something insensitive.</b>", 
    // name: "AMI_13", labels: respOptsAMI, required:true, horizontal: true},
    // {prompt: "<b>I start conversations without being prompted.</b>", 
    // name: "AMI_14", labels: respOptsAMI, required:true, horizontal: true},
    {prompt: "<b>When I have something I need to do, I do it straightaway so it is out of the way.</b>", 
    name: "AMI_15", labels: respOptsAMI, required:true, horizontal: true},
    // {prompt: "<b>I feel bad when I hear an acquaintance has an accident or illness.</b>", 
    // name: "AMI_16", labels: respOptsAMI, required:true, horizontal: true},
    // {prompt: "<b>I enjoy choosing what to do from a range of activities.</b>", 
    // name: "AMI_17", labels: respOptsAMI, required:true, horizontal: true},
    // {prompt: "<b>If I realise I have been unpleasant to someone, I will feel terribly guilty afterwards.</b>", 
    // name: "AMI_18", labels: respOptsAMI, required:true, horizontal: true}
  ],
  button_label: 'continue',
  scale_width: scaleDisplayWidth,
  on_finish: function(){ 
    var respData = jsPsych.data.getLastTrialData().select('response').values[0];
    var timeElapsed = jsPsych.getTotalTime();
    //console.log(respData); console.log(timeElapsed);  // for debugging only
    saveQuestData("AMI", respData, timeElapsed);
  }
};
// data_quest_wide$AMI_behavActiv <- rowMeans(AMI[c(5,9,10,11,12,15)]) # tendency to self-initiate goal-directed behaviour 
// data_quest_wide$AMI_socialMotiv <- rowMeans(AMI[c(2,3,4,8,14,17)])  # level of engagement in social interactions
// data_quest_wide$AMI_emoSens <- rowMeans(AMI[c(1,6,7,13,16,18)])     # feelings of positive and negative affection


////////////////////BADS-SF//////////////////////////////
// define response labels
var respOptsBADS = ["0    (Not at all)", "1", "2    (A little)", "3", "4       (A lot)", "5", "6 (Completely)"];

// scale items
var BADS = {
  type: jsPsychSurveyLikert,
  preamble: ("<p>"+
             "Please read each statement carefully and then select the number which best describes how "+
             "much the statement was true for you DURING THE PAST WEEK, INCLUDING TODAY."+
             "</p>"),
  questions: [
    {prompt: "<b>There were certain things I needed to do that I didn’t do.</b>", 
      name: "BADS_1", labels: respOptsBADS, required:true, horizontal: true}, 
    {prompt: "<b>I am content with the amount and types of things I did.</b>", 
      name: "BADS_2", labels: respOptsBADS, required:true, horizontal: true}, 
    {prompt: "<b>I engaged in many different activities.</b>", 
      name: "BADS_3", labels: respOptsBADS, required:true, horizontal: true},
    {prompt: "<b>I made good decisions about what type of activities and/or situations I put myself in.</b>", 
      name: "BADS_4", labels: respOptsBADS, required:true, horizontal: true}, 
    {prompt: "<b>I was an active person and accomplished the goals I set out to do.</b>", 
      name: "BADS_5", labels: respOptsBADS, required:true, horizontal: true}, 
    {prompt: "<b>Most of what I did was to escape from or avoid something unpleasant.</b>", 
      name: "BADS_6", labels: respOptsBADS, required:true, horizontal: true},
    {prompt: "<b>I spent a long time thinking over and over about my problems.</b>", 
      name: "BADS_7", labels: respOptsBADS, required:true, horizontal: true}, 
    {prompt: "<b>I engaged in activities that would distract me from feeling bad.</b>", 
      name: "BADS_8", labels: respOptsBADS, required:true, horizontal: true}, 
    {prompt: "<b>I did things that were enjoyable.</b>", 
      name: "BADS_9", labels: respOptsBADS, required:true, horizontal: true}
  ],
  button_label: 'continue',
  scale_width: scaleDisplayWidth,
  on_finish: function() {
    // // set progress bar manually
    // this.type.jsPsych.setProgressBar(0.8);
    // get response and RT data
    var respData = this.type.jsPsych.data.getLastTrialData().trials[0].response;
    var respRT = this.type.jsPsych.data.getLastTrialData().trials[0].rt;
    saveQuestData("BADS", respData, respRT);
  }
};


////////////////////ERQ-CR//////////////////////////////
// define response labels
var respOptsERQCR = ["1 (strongly disagree)", "2", "3", "4 (neutral)", "5", "6", "7 (strongly agree)"];

// scale items
var ERQCR = {
  type: jsPsychSurveyLikert,
  preamble: ("<p>"+
             "Finally, we would like to ask you some questions about your emotional life, in particular, how you "+
             "control (that is, regulate and manage) your emotions. Although some of the following "+
             "questions may seem similar to one another, they differ in important ways. For each statement "+
             "select how much you agree, thinking in particular about THE PAST WEEK."+
             "</p>"),
  questions: [
    {prompt: "<b>When I want to feel more <i>positive</i> emotion (such as joy or amusement), I change what I’m thinking about.</b>", 
      name: "ERQCR_1", labels: respOptsERQCR, required:true, horizontal: true}, 
    {prompt: "<b>When I want to feel less <i>negative</i> emotion (such as sadness or anger), I change what I’m thinking about.</b>", 
      name: "ERQCR_2", labels: respOptsERQCR, required:true, horizontal: true}, 
    {prompt: "<b>When I’m faced with a stressful situation, I make myself think about it in a way that helps me stay calm.</b>", 
      name: "ERQCR_3", labels: respOptsERQCR, required:true, horizontal: true},
    {prompt: "<b>When I want to feel more <i>positive</i> emotion, I change the way I’m thinking about the situation.</b>", 
      name: "ERQCR_4", labels: respOptsERQCR, required:true, horizontal: true},
    {prompt: "<b>I can speak over thirty languages fluently. </b>", 
      name: "catch_1", labels: respOptsERQCR, required:true, horizontal: true}, 
    {prompt: "<b>I control my emotions by changing the way I think about the situation I’m in.</b>", 
      name: "ERQCR_5", labels: respOptsERQCR, required:true, horizontal: true},
    {prompt: "<b>When I want to feel less <i>negative</i> emotion, I change the way I’m thinking about the situation.</b>", 
      name: "ERQCR_6", labels: respOptsERQCR, required:true, horizontal: true}
  ],
  button_label: 'continue',
  scale_width: scaleDisplayWidth,
  on_finish: function() {
    // // set progress bar manually
    // this.type.jsPsych.setProgressBar(0.8);
    // get response and RT data
    var respData = this.type.jsPsych.data.getLastTrialData().trials[0].response;
    var respRT = this.type.jsPsych.data.getLastTrialData().trials[0].rt;
    saveQuestData("ERQCR", respData, respRT);
  }
};



//////////////////////////// DEMOGS ///////////////////////
//age in years, gender identity, income bracket, employment status, housing status
var demogs = {
  type: jsPsychSurvey,
  pages: [
    [
      {
        type: 'html',
        prompt: "<b>Finally, we would like to ask you some questions about yourself "+
                "and your personal circumstances.</b>",
      },
      {
        type: 'text',
        prompt: "How old are you (in years)?", 
        name: 'demogs_age', 
        textbox_columns: 5,
        required: true,
        validation: "^[18-100]$"  // doesn't seem to work currently
      },
      {
        type: 'drop-down',
        prompt: "What is your gender identity?", 
        name: 'demogs_gender', 
        options: ['man', 'woman', 'non-binary', 'other', 'prefer not to say'], 
        required: true
      }, 
      {
        type: 'drop-down',
        prompt: "Which of the options below best describes your current employment status?", 
        name: 'demogs_employment', 
        options: ['employed (including full-time and part-time employment)', 
                  'unemployed (job seekers and those unemployed owing to ill health)',
                  'not seeking employment (stay-at-home parents, students, and retirees)'
                  ], 
        required: true
      },
      {
        type: 'drop-down',
        prompt: "Which of the options below best describes your current financial situation?", 
        name: 'demogs_financial', 
        options: ['doing okay financially',
                  'just about getting by',
                  'struggling financially'], 
        required: true
      },
      {
        type: 'drop-down',
        prompt: "Which of the options below best describes your current housing situation?", 
        name: 'demogs_housing', 
        options: ['homeowner (including those with a mortgage)',
                  'tenant',
                  'other (living with family or friends, homeless, or living in a hostel)'], 
        required: true
      },
      {
        type: 'multi-select',
        prompt: "Have you ever previously received treatment for a mental health problem? Please select all that apply",
        name: 'demogs_tx',
        options: ['yes - talking therapy (including cognitive-behavioural therapies)',
                  'yes - medication',
                  'yes - self-guided (e.g., workbooks or apps)',
                  'yes - other',
                  'no',
                  'prefer not to say'],
        required: true
      },
      {
        type: 'drop-down',
        prompt: "Do you consider yourself to be neurodivergent? "+
                "(Neurodivergence is a term for when someone processes or learns information in a different way to that which is considered 'typical': "+
                "common examples include autism and ADHD.)",
        name: 'demogs_neurodiv',
        options: ['yes', 'no', 'prefer not to say'],
        required: true
      },
      {
        type: 'multi-select',
        prompt: "Do you consider yourself to have a disability or form of neurodivergence that affects "+
                "your ability to do any of the below? Please select all that apply",
        name: 'demogs_disability',
        options: ['concentrate for extended periods of time',
                  'perform physically effortful activites',
                  'read, write, or do maths',
                  'deal with people you do not know',
                  'other form of impact not listed above',
                  'none of the above',
                  'prefer not to say'],
        required: true
      },
      {
        type: 'drop-down',
        prompt: "In the future, would you be willing to play more games like "+
                "the ones in this study, if you thought they could be used to "+
                "give you information about your thought processes or decision-making?", 
        options: ['yes', 'no', 'not sure'],
        //rows: 1,
        name: 'study_acceptability', 
        required: true
      }, 
      {
        type: 'text',
        prompt: "Over the last two weeks, did you experience any life events which "+
                "you think have significantly affected your answers to the questions "+
                "about your mood and feelings? For example, a change in your housing or "+
                "employment status, illness, or any other event? If you want to, you "+
                "can tell us about this below (optional question).", 
        name: 'life_events', 
        textbox_rows: 4,
        textbox_columns: 60,
        required: false
      },
      {
        type: 'text',
        prompt: "Is there any feedback you would like to give us about "+
                "any aspect of the study (including the tasks and questionnaires)?", 
        name: 'study_feedback', 
        textbox_rows: 4,
        textbox_columns: 60,
        required: false
      },
    ]
  ],
  show_question_numbers: 'onPage',
  button_label_finish: 'submit',
  on_finish: function() {
    // get response and RT data
    var respData = this.type.jsPsych.data.getLastTrialData().trials[0].response;
    var respRT = this.type.jsPsych.data.getLastTrialData().trials[0].rt;
    saveQuestData("demogs", respData, respRT);
    // end of study
    saveEndData();
  }
};

///////////////////////////////////////////// CONCAT /////////////////////////////////////////////////////////
var timeline_quests = [];      
timeline_quests.push(questsIntroText);
timeline_quests.push(PHQ9);
timeline_quests.push(BADS);
timeline_quests.push(AMI);
timeline_quests.push(DAS);
timeline_quests.push(miniSPIN);
timeline_quests.push(ERQCR);
timeline_quests.push(demogs);
timeline_quests.push(questsEndScreen);

export { timeline_quests };

///////////////////////////////////////////// EXPORT MODULAR ///////////////////////////////////////////////////
// function to run self-reports separately for experiments coming from phaser
export function runQuests() {
  // run timelines
  jsPsych.run(timeline_quests); 
};

