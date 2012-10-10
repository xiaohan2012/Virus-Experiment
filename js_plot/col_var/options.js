var seriesDefaultsOptions = {
    renderer:$.jqplot.BarRenderer,
    rendererOptions: {
        barWidth: 50,
    }
};
var step= 0.1;
var ticks = [];
/*var seriesOptions = [];*/
for(var i= 0;i<= 1; i+=step){
    ticks.push(i.toFixed(2).toString());
    /*seriesOptions .push({label:i.toFixed(2).toString()});*/
}

var axesOptions = {
            xaxis:{
                min: 0,
                max: 1.0,
                label:'Similarity score group',
                labelOptions:{
                    fontSize:18
                },
                ticks:ticks,
                tickOptions:{
                    fontSize:18,
                }
            },
            yaxis:{
                min: 0.0,
                label:'Frequency',
                labelRenderer: $.jqplot.CanvasAxisLabelRenderer,
                tickOptions:{
                    fontSize:18,
                }
            } 
};
