var maxSize = 0;
var maxCoverage = 0.0;
var mapMargin = 30;
var mismatchRatioThreshold = 0.05;
var sorting_by_coverage_rate = 1;

function switch_sort() {
    if(sorting_by_coverage_rate == 1) {
        sorting_by_coverage_rate = 0;
        genome_coverage_data.sort(sortByBases);
        drawCoverages('genome_coverage', genome_coverage_data, genome_sizes, stats_bin);
        document.getElementById("sort_by_div").innerHTML = "Order by: <a href='javascript:switch_sort();'>Coverage rate</a> | <font color='#FF6600'>Bases on target</font>";
    } else {
        sorting_by_coverage_rate = 1;
        genome_coverage_data.sort(sortByCoverageRate);
        drawCoverages('genome_coverage', genome_coverage_data, genome_sizes, stats_bin);
        document.getElementById("sort_by_div").innerHTML = "Order by: <font color='#FF6600'>Coverage rate</font> | <a href='javascript:switch_sort();'>Bases on target</a>";
    }
}

function sortByCoverageRate(a, b) {
    var r = b["coverage_rate"] -  a["coverage_rate"];
    if(r == 0)
        return b["bases"] -  a["bases"];
    else
        return r;
}

function sortByBases(a, b) {
    var r= b["bases"] -  a["bases"];
    if(r == 0)
        return b["coverage_rate"] -  a["coverage_rate"];
    else
        return r;
}

function drawCoverages(divid, data, sizes, bin) {
    mapcontainer = document.getElementById(divid);
    var hasData = false;
    for(d in data) {
        hasData = true;
        if(sizes[d] > maxSize)
            maxSize = sizes[d];

        for(c in data[d]["coverage"]) {
            if(data[d]["coverage"][c] > maxCoverage)
                maxCoverage = data[d]["coverage"][c];
        }
    }
    if(hasData) {
        mapcontainer.style.display = 'block';
    } else {
        mapcontainer.style.display = 'none';
    }
    var childs = mapcontainer.childNodes;
    for(var i = 0; i < childs.length; i++) {
      mapcontainer.removeChild(childs[i]);
    }
    var colorTableHTML = "<div id='colortable' style='text-align:center;padding:5px;font-size:14px;font-family:Arial;'>";
    colorTableHTML += "<center><table style='border:0px;'> <tr>  <td style='width:200px; border:0px;text-align:right;'> Mismatch ratio = 0 </td>";
    var count = 0;
    var tds = 30;
    while (count < tds) {
        var mr = mismatchRatioThreshold * count/tds ;
        count++;
        var c = getColor(mr);
        colorTableHTML += "<td style='background:" + c + "; width=10px;' title = '" + mr + "'></td>  "; 
    }
    colorTableHTML += "<td style='width:200px;border:0px;'> Mismatch ratio = " + mismatchRatioThreshold + " </td>";
    colorTableHTML += "</tr></table></center></div>";
    mapcontainer.innerHTML = colorTableHTML;

    for(d in data) {
        var genome = data[d];
        cvs = document.createElement("canvas");
        cvs.id = 'coverage_' + d.toString();
        cvs.width=mapcontainer.offsetWidth - 10;
        cvs.height=60;
        cvs.style.padding='2px 0px 2px 0px';
        cvs.onmousemove = onCanvasMove;
        cvs.onmouseover = onCanvasIn;
        cvs.onmouseout = onCanvasOut;
        mapcontainer.appendChild(cvs);

        drawGenome(genome, cvs.id, sizes[d], bin);

        var namediv = document.createElement("div"); 
        namediv.innerHTML = "<div style='text-align:center;padding:2px;font-size:10px;color:#339967;'> " + data[d]['name'] + " (" + data[d]['coverage_rate'] + "% covered, " + data[d]['reads'] + " reads, " + data[d]['bases'] + " bases)</div>" ;; 
        mapcontainer.appendChild(namediv); 
    }
}

function onCanvasMove(e) {
    var cvs = e.target;
    var id = parseInt(cvs.id.substring(9));
    console.log(cvs.id.substring(9));
    console.log(id);
    if(!genome_coverage_data[id])
        return;

    genome = genome_coverage_data[id];
    var mapw = cvs.width;
    var maph = cvs.height;
    var x = e.clientX - cvs.offsetLeft;
    var pos = (x - mapMargin) * maxSize / (mapw - mapMargin*2);
    pos = Math.round(pos / stats_bin);

    console.log(pos);

    if(!genome["coverage"][pos])
        return;

    var start = pos * stats_bin;
    var end = (pos+1) * stats_bin - 1;
    var html = start + "-" + end + "<br>";
    html += "mean coverage: " + genome["coverage"][pos] + "<br>";
    html += "mean mismatch ratio: " + genome["mismatch_ratios"][pos] + "<br>";

    var tips = document.getElementById('maptips');
    tips.style.position="absolute";
    tips.style.left = e.clientX + 5 + tips.parentNode.scrollLeft;
    tips.style.top = e.clientY + 5 + tips.parentNode.scrollTop + document.body.scrollTop;
    tips.innerHTML = html;

    console.log(html);
}

function onCanvasIn() {
    var tips = document.getElementById('maptips');
    tips.style.display = 'block';
}

function onCanvasOut() {
    var tips = document.getElementById('maptips');
    tips.style.display = 'none';
}

function drawGenome(genome, canvasid, size, bin) {
    var cvs = document.getElementById(canvasid);
    var mapw = cvs.width;
    var maph = cvs.height;
    var ctx = cvs.getContext("2d"); 
    var texth = 15;

    var name = genome['name'];
    var reads = genome['reads'];
    var bases = genome['bases'];
    var avg_mismatch_ratio = genome['avg_mismatch_ratio'];
    var coverage = genome['coverage'];
    var mismatch_ratios = genome['mismatch_ratios'];

    for(var pos in coverage) {
        var c = coverage[pos];
        var mr = mismatch_ratios[pos];

        var w = (mapw - mapMargin*2) * bin / maxSize;
        var drawW = w;
        var h = (maph-texth)* (c/maxCoverage);
        var drawH = h;
        var centerX = mapMargin + (pos-0.5) * (mapw - mapMargin*2) * bin/maxSize;
        var x = centerX - 1;
        var y = maph- texth;
        ctx.fillStyle=getColor(mr);
        ctx.fillRect(x, y, drawW, -drawH);
    }

    var xbars = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    var tailPainted = false;
    for(b in xbars) {
        var tick = Math.round( maxSize * xbars[b]/10 );
        if(tick > size) {
            if(tailPainted)
                continue;
            tailPainted = true;
            tick = size;
        }
        var x = mapMargin + (mapw - mapMargin*2) * tick/maxSize;
        x = Math.round(x);
        ctx.font = "10px Arial";
        ctx.fillStyle = "#999999";
        ctx.fillText(tick.toString(), x, maph);
    }

    ctx.font = "10px Arial";
    ctx.fillStyle = "#AAAAAA";
    ctx.fillText(maxCoverage.toString() + "", 10, 10);
}

function getColor(mr) {
    if(mr > mismatchRatioThreshold)
        return "rgb(128, 128, 128)";
    else {
        var p = mr/mismatchRatioThreshold;
        var diff = 120* p;
        return "rgb(" + (128+diff) + "," + (128-diff) + "," + (128-diff) + ")";
    }
}