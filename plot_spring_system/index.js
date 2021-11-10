import { produceSystemSvg } from './produce_system_svg.js';

const [ filename1, filename2, outputFilename ] = process.argv.slice(2);

produceSystemSvg(filename1, filename2, outputFilename).then(() => console.log("Done"));
