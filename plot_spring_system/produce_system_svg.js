import { JSDOM } from "jsdom";

import * as d3 from "d3";
import * as fs from "fs";

const SVG_WIDTH = 800;
const SVG_HEIGHT = 800;

let svg;

function initialize(nodes, links) {
    // Initialize the nodes
    var link = svg
        .selectall("line")
        .data(links)
        .enter()
        .append("line")
        .style("stroke", "#aaa")
        .style("stroke-width", (d, i) => d.value)
        .attr("x1", (d, i) => nodes[d.source].x)
        .attr("y1", (d, i) => nodes[d.source].y)
        .attr("x2", (d, i) => nodes[d.target].x)
        .attr("y2", (d, i) => nodes[d.target].y);
    var node = svg
        .selectAll("circle")
        .data(nodes)
        .enter()
        .append("circle")
        .attr("r", 20)
        .attr("cx", (d, i) => d.x)
        .attr("cy", (d, i) => d.y)
        .style("stroke", "#000")
        .style("fill", "#69b3a2");
}

function animate(nodes, links) {
    // Initialize the nodes
    var duration = 1000;
    var link = svg
        .selectAll("line")
        .data(links)
        .transition()
        .duration(duration)
        .style("stroke", "#aaa")
        .style("stroke-width", (d, i) => d.value)
        .attr("x1", (d, i) => nodes[d.source].x)
        .attr("y1", (d, i) => nodes[d.source].y)
        .attr("x2", (d, i) => nodes[d.target].x)
        .attr("y2", (d, i) => nodes[d.target].y);
    svg.selectAll("circle")
        .data(nodes)
        .transition()
        .duration(duration)
        .attr("cx", (d, i) => d.x)
        .attr("cy", (d, i) => d.y);
}

var g_nodes = null;
var g_links = null;

const config = {
    radius: 10,
    padding: 20,
};
function setData(nodes, links) {
    // remove old content:
    // d3.selectAll("svg > *").remove();

    console.log("number of nodes: ", nodes.length);
    console.log("number of links: ", links.length);
    if (nodes[0].length != 2) {
        throw new Error("node has wrong number of columns, should be: x,y");
    }
    if (links[0].length != 3) {
        throw new Error(
            "link has wrong number of columns, should be: src,dst,weight"
        );
    }
    links.forEach(function (link) {
        link[0] = parseInt(link[0]);
        link[1] = parseInt(link[1]);
    });
    var min_x = Number.MAX_VALUE;
    var min_y = Number.MAX_VALUE;
    var max_x = Number.MIN_VALUE;
    var max_y = Number.MIN_VALUE;
    nodes.forEach(function (node) {
        node[0] = parseFloat(node[0]);
        node[1] = parseFloat(node[1]);
        min_x = Math.min(min_x, node[0]);
        min_y = Math.min(min_y, node[1]);
        max_x = Math.max(max_x, node[0]);
        max_y = Math.max(max_y, node[1]);
    });
    console.log("x range", min_x, max_x);
    console.log("y range", min_y, max_y);
    const scaleX = d3
        .scaleLinear()
        .domain([min_x, max_x])
        .range([config.padding, SVG_WIDTH - config.padding]);
    const scaleY = d3
        .scaleLinear()
        .domain([min_y, max_y])
        .range([config.padding, SVG_HEIGHT - config.padding]);

    nodes.forEach(function (node) {
        node[0] = scaleX(node[0]);
        (node[0] < 0 || node[0] > SVG_WIDTH) && console.error("scaled out of viewport:", node[0]);
        node[1] = scaleY(node[1]);
        (node[1] < 0 || node[1] > SVG_HEIGHT) && console.error("scaled out of viewport:", node[1]);
    });
    // console.log(nodes);

    g_nodes = nodes;
    g_links = links;

    var link = svg
        .selectAll("line")
        .data(links)
        .enter()
        .append("line")
        .style("stroke", "#aaa")
        .style("stroke-width", (d, i) => d[2])
        .attr("x1", (d, i) => nodes[d[0]][0])
        .attr("y1", (d, i) => nodes[d[0]][1])
        .attr("x2", (d, i) => nodes[d[1]][0])
        .attr("y2", (d, i) => nodes[d[1]][1]);
    var node = svg
        .selectAll("circle")
        .data(nodes)
        .enter()
        .append("circle")
        .attr("r", config.radius)
        .attr("cx", (d, i) => d[0])
        .attr("cy", (d, i) => d[1])
        .style("stroke", "#000")
        .style("fill", "#69b3a2");
}

function readFileCSV(filename, callback) {
    const content = fs.readFileSync(filename).toString();
    // const reader = new FileReader();
    // reader.readAsText(filename);
    // reader.onloadend = function (ev) {
    //     if (reader.readyState !== FileReader.DONE) {
    //         console.log("FileReader failed for " + filename.name);
    //         return;
    //     }
    const data = d3.csvParseRows(content);
    callback(data);
    // };
}

function isNodesFilename(fn) {
    return fn.search(/nodes/i) !== -1;
}

function isLinksFilename(fn) {
    return fn.search(/links|edges/i) !== -1;
}

function filenameType(fn) {
    const n = isNodesFilename(fn);
    const l = isLinksFilename(fn);
    if (!n && !l) {
        throw new Error(
            "could not identify filename '" + fn + "' as nodes or edges file"
        );
    }
    if (n && l) {
        throw new Error("ambiguous filename '" + fn + "'");
    }
    if (n) return "nodes";
    if (l) return "links";
}
function handleFiles(filenames) {
    if (filenames.length != 2) {
        throw new Error("Select two files - nodes.csv and links.csv");
    }

    const t0 = filenameType(filenames[0]);
    const t1 = filenameType(filenames[1]);
    if (t0 === null || t1 === null) return;
    if (t0 === t1) {
        throw new Error("select one(!) node file and one(!) links file");
    }
    if (t0 == "nodes") {
        var nodesFile = filenames[0];
        var linksFile = filenames[1];
    } else {
        var nodesFile = filenames[1];
        var linksFile = filenames[0];
    }
    readFileCSV(nodesFile, function (nd) {
        readFileCSV(linksFile, function (ld) {
            setData(nd, ld);
        });
    });
}

/**
 * Mostly copied from https://slides.cgg.unibe.ch/aopt20/plots/springsys.html
 *
 * @param {string} filename1 path + filename
 * @param {string} filename2 path + filename
 * @param {string} svgOutput path + filename
 */
export function produceSystemSvg(filename1, filename2, svgOutput) {
    // svg.append("rect")
    //     .attr("x", 10)
    //     .attr("y", 10)
    //     .attr("width", 80)
    //     .attr("height", 80)
    //     .style("fill", "orange");
    // svg.append("rect")
    //     .attr("x", 10)
    //     .attr("y", 10)
    //     .attr("width", 20)
    //     .attr("height", 20)
    //     .style("fill", "red");

    const dom = new JSDOM(`<!DOCTYPE html><body style="margin: 0"></body>`);

    const body = d3.select(dom.window.document.querySelector("body"));

    svg = body
        .append("svg")
        .attr("width", SVG_WIDTH)
        .attr("height", SVG_HEIGHT)
        .attr("xmlns", "http://www.w3.org/2000/svg");

    handleFiles([filename1, filename2]);
    // setupDropzone(document, handleFiles);

    fs.writeFileSync(svgOutput, body.html());
}
