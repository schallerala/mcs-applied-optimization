import { join, basename } from 'path';
import { globbySync } from 'globby';

import { produceSystemSvg } from './produce_system_svg.js';

const [ resultOutput ] = process.argv.slice(2);

const allCsv = globbySync(join(resultOutput, '*.csv')).sort();
// as sorted, can expect the results to be ordered with
// - 2 csv for before
// - 2 csv for after (optimized version)
if (allCsv.length % 4 != 0)
    throw new Error("Expected a multiple of 4 csv files (2 before 2 after: nodes and links) - Got: " + allCsv.length);

(async function () {
    for (let i = 0; i < allCsv.length; i += 4) {
        const file1Before = allCsv[i];
        const file2Before = allCsv[i + 1];

        const file1After = allCsv[i + 2];
        const file2After = allCsv[i + 3];

        // ex: <path>/execution_GradientDescent_0_0_10__1_edges.csv
        const outputFilenameBase = file1Before.slice(0, -"_edges.csv".length);

        await produceSystemSvg(file1Before, file2Before, `${outputFilenameBase}_before`);
        await produceSystemSvg(file1After, file2After, `${outputFilenameBase}_optimized`);
    }
})().then(() => console.log("Done"));