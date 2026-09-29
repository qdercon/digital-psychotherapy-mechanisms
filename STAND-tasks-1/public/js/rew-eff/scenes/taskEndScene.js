// End scene to inform participants they have finished the task, and route them to the post-task questions

// import task info from versionInfo file
import { taskCond } from "../../versionInfo.js"; 

// import js game element modules (sprites, ui, outcome animations)
import InstructionsPanel from "../elements/instructionsPanel.js";

// import our custom events centre for passsing info between scenes and data saving function
import eventsCenter from "../eventsCenter.js";

// import functions to run next activity
import { runStudyCaus } from "../../causal-attr/constructStudy.js";
import { runQuests } from "../../selfReports.js";

// this function extends Phaser.Scene and includes the core logic for the scene
export default class TaskEndScene extends Phaser.Scene {
    constructor() {
        super({
            key: 'TaskEndScene'
        });
    }

    preload() {
        // load cloud sprites to add texture to background
        this.load.image('cloud1', '../../assets/imgs/cloud1.png');
    }
    
    create() {
        // load a few cloud sprites dotted around
        const cloud1 = this.add.sprite(180, 100, 'cloud1');
        const cloud2 = this.add.sprite(320, 540, 'cloud1');
        const cloud3 = this.add.sprite(630, 80, 'cloud1');
        var gameHeight = this.sys.game.config.height;
        var gameWidth = this.sys.game.config.width;

        var titleText = 'Game Over!'
        ///////////////////PAGE ONE////////////////////
        var mainTxt = ("Thank you for playing.\n\n" +
                        "Press the button below to continue to the next activity.\n\n");
        var buttonTxt = "continue";
        var pageNo = 1;
        this.endPanel = new InstructionsPanel(this, gameWidth/2, gameHeight/2,
                                                pageNo, titleText, mainTxt, buttonTxt);
        // end scene
        eventsCenter.once('page1complete', function () {
            // whatever condition, we'll be moving back to jspsych next
            document.getElementById('game-container').style.display = "none";   // hide phaser container
            document.getElementById('jspsych-target').style.display = "block";  // show jspsych container
            document.body.style.background = "aliceblue";
            // do next part
            this.nextScene();
        }, this);
    }
    
    update(time, delta) {
    }
    
    nextScene() {
        // if this was the first task, move to the second task, else go to self-reports
        if ( taskCond == "rew-eff-first") {
            runStudyCaus();
        } else {
            runQuests();
        }
    } 
}