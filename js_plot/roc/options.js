var seriesDefaultsOptions = {lineWidth: 2,showMarker:false};
var legendOptions = {show:true,location:"se",labels:["1","2","3","4","5","6"],fontSize:"18px"};
var axesOptions = {
            xaxis:{
                labelRenderer: $.jqplot.CanvasAxisLabelRenderer,
                min: 0.0,
                max: 1.0,
                tickInterval:0.1,
                tickOptions:{
                    formatString:"%1.1f",
                    fontSize: "16px",
                },
                label:'1-specificity',
                labelOptions:{
                    fontSize:"22px",
                }
            }, 
            yaxis:{
                min: 0,
                max: 1.0,
                tickOptions:{
                    formatString:"%1.1f",
                    fontSize: "16px",
                },
                label:'sentivity',
                labelRenderer: $.jqplot.CanvasAxisLabelRenderer,
                labelOptions:{
                    fontSize:"22px",
                }
            }
};

