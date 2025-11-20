/**
 * Extract leaderboard table from PolarisHub - Run in Browser Console
 *
 * Instructions:
 * 1. Navigate to: https://polarishub.io/competitions/asap-discovery/antiviral-drug-discovery-2025
 * 2. Open Developer Tools (F12)
 * 3. Go to Console tab
 * 4. Copy and paste this entire script
 * 5. Press Enter
 * 6. Data will be copied to clipboard and logged to console
 */

(function extractLeaderboard() {
    console.log('🔍 Searching for leaderboard table...');

    // Try to find the table
    let table = document.querySelector('table');

    if (!table) {
        console.error('❌ No table found on page. Make sure the leaderboard is visible.');
        return;
    }

    // Extract headers
    let headers = Array.from(table.querySelectorAll('th')).map(th => th.innerText.trim());

    if (headers.length === 0) {
        // Try alternative: first row might be headers
        let firstRow = table.querySelector('tr');
        headers = Array.from(firstRow.querySelectorAll('td, th')).map(cell => cell.innerText.trim());
    }

    console.log('📊 Headers found:', headers);

    // Extract rows
    let rows = Array.from(table.querySelectorAll('tr')).slice(1); // Skip header row

    let data = rows.map(row => {
        let cells = Array.from(row.querySelectorAll('td'));
        let rowData = {};

        cells.forEach((cell, index) => {
            let headerName = headers[index] || `Column_${index}`;
            rowData[headerName] = cell.innerText.trim();
        });

        return rowData;
    }).filter(row => Object.keys(row).length > 0); // Remove empty rows

    console.log(`✅ Extracted ${data.length} rows`);
    console.log('Preview:', data.slice(0, 3));

    // Convert to CSV format
    let csv = headers.join(',') + '\n';
    data.forEach(row => {
        let values = headers.map(h => {
            let val = row[h] || '';
            // Escape commas and quotes
            if (val.includes(',') || val.includes('"')) {
                val = '"' + val.replace(/"/g, '""') + '"';
            }
            return val;
        });
        csv += values.join(',') + '\n';
    });

    // Copy JSON to clipboard
    let jsonData = JSON.stringify(data, null, 2);

    // Try to copy to clipboard
    if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(jsonData).then(() => {
            console.log('📋 JSON data copied to clipboard!');
        }).catch(err => {
            console.log('⚠️  Could not copy to clipboard automatically');
        });
    }

    // Also copy CSV to a temporary textarea
    let csvTextarea = document.createElement('textarea');
    csvTextarea.value = csv;
    document.body.appendChild(csvTextarea);
    csvTextarea.select();
    try {
        document.execCommand('copy');
        console.log('📋 CSV data also available (copied)');
    } catch (err) {
        console.log('⚠️  Could not copy CSV');
    }
    document.body.removeChild(csvTextarea);

    // Log results
    console.log('\n📄 JSON Format:');
    console.log(jsonData);
    console.log('\n📄 CSV Format:');
    console.log(csv);

    // Download as file
    function downloadFile(content, filename, contentType) {
        let blob = new Blob([content], { type: contentType });
        let url = URL.createObjectURL(blob);
        let link = document.createElement('a');
        link.href = url;
        link.download = filename;
        link.click();
        URL.revokeObjectURL(url);
    }

    // Offer downloads
    console.log('\n💾 To download files, run:');
    console.log('downloadAsJSON() - Download as JSON');
    console.log('downloadAsCSV() - Download as CSV');

    window.downloadAsJSON = function() {
        downloadFile(jsonData, 'leaderboard.json', 'application/json');
        console.log('✅ Downloaded leaderboard.json');
    };

    window.downloadAsCSV = function() {
        downloadFile(csv, 'leaderboard.csv', 'text/csv');
        console.log('✅ Downloaded leaderboard.csv');
    };

    return {
        data: data,
        csv: csv,
        json: jsonData,
        downloadJSON: window.downloadAsJSON,
        downloadCSV: window.downloadAsCSV
    };
})();
