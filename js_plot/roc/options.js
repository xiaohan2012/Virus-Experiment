var seriesDefaultsOptions = {lineWidth: 2,showMarker:false};
var legendOptions = {show:true,location:"se",labels:["1","2","3","4","5","6"],fontSize:"18px"};
var axesOptions = {
            xaxis:{
                min: 0.0,
                max: 1.0,
                tickOptions:{
                    formatString:"%1.1f",
                },
                label:'specificity',
                labelRenderer: $.jqplot.CanvasAxisLabelRenderer,
            }, 
            yaxis:{
                min: 0,
                max: 1.0,
                tickOptions:{
                    formatString:"%1.1f",
                },
                label:'sentivity',
                labelRenderer: $.jqplot.CanvasAxisLabelRenderer,
            }
};

