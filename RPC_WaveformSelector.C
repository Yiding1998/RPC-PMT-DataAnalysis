// RPC_WaveformSelector.C
//
// Run in ROOT, for example:
//   root -l -b -q 'RPC_WaveformSelector.C(20,40,0.3,0.8,0.8,1.5)'
//
// The CSV format is three columns: time, rpc, pmt.
// Time is converted from s to ns. RPC/PMT signals are converted from V to mV
// and inverted so negative pulses become positive.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TTree.h"

namespace fs = std::filesystem;
using namespace std;

const string kDefaultCsvFolder =
    "/ustcfs/STCFUser/yzhao/LowPressure/ExpDatas_2606/40_50_10/3kpa";
const double kBaselineWindowNs = 10.0;

struct RpcShapeParams {
    double amplitude;
    double rise_time;
    double fwhm;
    int amplitude_idx;
    bool valid;
};

double interpolateCrossingTime(const vector<double>& time,
                               const vector<double>& signal,
                               double threshold,
                               int idx1,
                               int idx2) {
    if (idx1 < 0 || idx2 < 0 ||
        idx1 >= static_cast<int>(signal.size()) ||
        idx2 >= static_cast<int>(signal.size())) {
        return 0.0;
    }

    const double y1 = std::abs(signal[idx1]);
    const double y2 = std::abs(signal[idx2]);
    const double t1 = time[idx1];
    const double t2 = time[idx2];

    if (std::abs(y2 - y1) < 1e-15) {
        return 0.5 * (t1 + t2);
    }

    return t1 + (threshold - y1) * (t2 - t1) / (y2 - y1);
}

vector<int> getVoltageFolders(const string& dirPath) {
    vector<int> voltages;

    TSystemDirectory dir(dirPath.c_str(), dirPath.c_str());
    TList* files = dir.GetListOfFiles();
    if (!files) {
        cerr << "Cannot open directory: " << dirPath << endl;
        return voltages;
    }

    TIter next(files);
    TSystemFile* file = nullptr;
    while ((file = static_cast<TSystemFile*>(next()))) {
        TString name = file->GetName();
        if (name == "." || name == "..") {
            continue;
        }
        if (file->IsDirectory() && name.IsDigit()) {
            voltages.push_back(name.Atoi());
        }
    }

    sort(voltages.begin(), voltages.end());
    delete files;
    return voltages;
}

int countCSVFiles(const string& folderPath) {
    int count = 0;
    try {
        for (const auto& entry : fs::directory_iterator(folderPath)) {
            if (entry.is_regular_file() && entry.path().extension() == ".csv") {
                ++count;
            }
        }
    } catch (const fs::filesystem_error& err) {
        cerr << "Cannot access folder " << folderPath << ": " << err.what() << endl;
        return 0;
    }
    return count;
}

bool readCSVData(const string& filename,
                 vector<double>& time,
                 vector<double>& rpc_signal,
                 vector<double>& pmt_signal) {
    ifstream file(filename, ios::binary);
    if (!file.is_open()) {
        return false;
    }

    time.clear();
    rpc_signal.clear();
    pmt_signal.clear();
    time.reserve(7000);
    rpc_signal.reserve(7000);
    pmt_signal.reserve(7000);

    string line;
    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        const char* ptr = line.c_str();
        char* end_ptr = nullptr;

        double t = strtod(ptr, &end_ptr) * 1e9;
        if (end_ptr == ptr) {
            continue;
        }

        ptr = end_ptr + 1;
        double rpc = -strtod(ptr, &end_ptr) * 1e3;
        if (end_ptr == ptr) {
            continue;
        }

        ptr = end_ptr + 1;
        double pmt = -strtod(ptr, &end_ptr) * 1e3;
        if (end_ptr == ptr) {
            continue;
        }

        time.push_back(t);
        rpc_signal.push_back(rpc);
        pmt_signal.push_back(pmt);
    }

    return !time.empty() && time.size() == rpc_signal.size() && time.size() == pmt_signal.size();
}

void baselineCorrection(vector<double>& signal,
                        const vector<double>& time,
                        double baseline_window_ns = kBaselineWindowNs) {
    if (signal.size() < 3 || time.size() < 3) {
        return;
    }

    const double step = time[2] - time[1];
    if (step <= 0) {
        return;
    }

    const double baseline_end = time[0] + baseline_window_ns;
    int end_idx = static_cast<int>((baseline_end - time[0]) / step);
    end_idx = std::max(1, std::min(end_idx, static_cast<int>(signal.size())));

    double sum = 0.0;
    for (int i = 0; i < end_idx; ++i) {
        sum += signal[i];
    }
    const double baseline = sum / end_idx;

    for (double& value : signal) {
        value -= baseline;
    }
}

RpcShapeParams calculateRpcShapeParams(const vector<double>& time,
                                       const vector<double>& signal) {
    RpcShapeParams params = {0.0, 0.0, 0.0, -1, false};
    if (time.size() < 3 || signal.size() != time.size()) {
        return params;
    }

    const double time_min = std::max(-14.0, time.front());
    const double time_max = std::min(14.0, time.back());
    const double step = time[2] - time[1];
    if (step <= 0) {
        return params;
    }

    int min_idx = static_cast<int>((time_min - time.front()) / step);
    int max_idx = static_cast<int>((time_max - time.front()) / step);
    min_idx = std::max(0, min_idx);
    max_idx = std::min(static_cast<int>(signal.size()) - 1, max_idx);
    if (max_idx <= min_idx) {
        return params;
    }

    auto max_it = max_element(signal.begin() + min_idx, signal.begin() + max_idx + 1);
    auto min_it = min_element(signal.begin() + min_idx, signal.begin() + max_idx + 1);
    const double max_val = *max_it;
    const double min_val = *min_it;
    const int max_amp_idx = static_cast<int>(max_it - signal.begin());
    const int min_amp_idx = static_cast<int>(min_it - signal.begin());

    params.amplitude = std::max(std::abs(max_val), std::abs(min_val));
    params.amplitude_idx = (std::abs(max_val) >= std::abs(min_val)) ? max_amp_idx : min_amp_idx;
    if (params.amplitude <= 0.0) {
        return params;
    }

    const double threshold_10 = 0.1 * params.amplitude;
    const double threshold_90 = 0.9 * params.amplitude;
    const double threshold_50 = 0.5 * params.amplitude;

    int idx_10_left = -1;
    int idx_90_left = -1;
    int idx_50_left = -1;
    int idx_50_right = -1;

    for (int i = params.amplitude_idx; i > 0; --i) {
        const double value = std::abs(signal[i]);
        if (idx_90_left == -1 && value <= threshold_90) {
            idx_90_left = i;
        }
        if (idx_50_left == -1 && value <= threshold_50) {
            idx_50_left = i;
        }
        if (idx_10_left == -1 && value <= threshold_10) {
            idx_10_left = i;
        }
        if (idx_90_left != -1 && idx_50_left != -1 && idx_10_left != -1) {
            break;
        }
    }

    for (int i = params.amplitude_idx; i < static_cast<int>(signal.size()); ++i) {
        if (std::abs(signal[i]) <= threshold_50) {
            idx_50_right = i;
            break;
        }
    }

    if (idx_10_left > 0 && idx_90_left > 0) {
        const double t10 = interpolateCrossingTime(time, signal, threshold_10,
                                                   idx_10_left - 1, idx_10_left);
        const double t90 = interpolateCrossingTime(time, signal, threshold_90,
                                                   idx_90_left - 1, idx_90_left);
        params.rise_time = std::abs(t90 - t10);
    }

    if (idx_50_left > 0 && idx_50_right > 0 && idx_50_right + 1 < static_cast<int>(signal.size())) {
        const double t_left = interpolateCrossingTime(time, signal, threshold_50,
                                                      idx_50_left - 1, idx_50_left);
        const double t_right = interpolateCrossingTime(time, signal, threshold_50,
                                                       idx_50_right + 1, idx_50_right);
        params.fwhm = std::abs(t_right - t_left);
    }

    params.valid = params.rise_time > 0.0 && params.fwhm > 0.0;
    return params;
}

bool inRange(double value, double min_value, double max_value) {
    return value >= min_value && value <= max_value;
}

TCanvas* makeWaveformCanvas(int voltage,
                            const string& csv_stem,
                            const vector<double>& time,
                            const vector<double>& rpc_signal,
                            const vector<double>& pmt_signal,
                            const RpcShapeParams& params) {
    TCanvas* canvas = new TCanvas(Form("c_%dV_%s", voltage, csv_stem.c_str()),
                                  Form("%dV %s", voltage, csv_stem.c_str()),
                                  1200, 800);
    canvas->SetGrid();

    TGraph* gr_rpc = new TGraph(static_cast<int>(time.size()));
    TGraph* gr_pmt = new TGraph(static_cast<int>(time.size()));
    gr_rpc->SetName(Form("g_rpc_%dV_%s", voltage, csv_stem.c_str()));
    gr_pmt->SetName(Form("g_pmt_%dV_%s", voltage, csv_stem.c_str()));
    gr_rpc->SetTitle("RPC");
    gr_pmt->SetTitle("PMT");

    for (int i = 0; i < static_cast<int>(time.size()); ++i) {
        gr_rpc->SetPoint(i, time[i], rpc_signal[i]);
        gr_pmt->SetPoint(i, time[i], pmt_signal[i]);
    }

    gr_rpc->SetLineColor(kBlue + 1);
    gr_rpc->SetLineWidth(2);
    gr_pmt->SetLineColor(kRed + 1);
    gr_pmt->SetLineWidth(2);

    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle(Form("%dV %s;Time (ns);Signal (mV)", voltage, csv_stem.c_str()));
    mg->Add(gr_rpc, "L");
    mg->Add(gr_pmt, "L");
    mg->Draw("A");

    TLegend* legend = new TLegend(0.42, 0.72, 0.90, 0.90);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.038);
    legend->SetMargin(0.12);
    legend->AddEntry(gr_rpc,
                     Form("RPC: A=%.3f mV, RT=%.3f ns, FWHM=%.3f ns",
                          params.amplitude, params.rise_time, params.fwhm),
                     "l");
    legend->AddEntry(gr_pmt, "PMT", "l");
    legend->Draw();

    canvas->Modified();
    canvas->Update();
    return canvas;
}

void RPC_WaveformSelector(double amp_min,
                          double amp_max,
                          double rise_min,
                          double rise_max,
                          double fwhm_min,
                          double fwhm_max,
                          string csvFolder = kDefaultCsvFolder,
                          const char* outputRootName = "selected_waveforms.root",
                          const char* pngBaseDir = "selected_waveforms_png",
                          int maxFilesPerVoltage = -1) {
    if (amp_min > amp_max || rise_min > rise_max || fwhm_min > fwhm_max) {
        cerr << "Invalid range: min value is larger than max value." << endl;
        return;
    }

    vector<int> voltages = getVoltageFolders(csvFolder);
    if (voltages.empty()) {
        cerr << "No numeric voltage folders found under " << csvFolder << endl;
        return;
    }

    fs::create_directories(pngBaseDir);

    TFile* output = new TFile(outputRootName, "RECREATE");
    if (!output || output->IsZombie()) {
        cerr << "Cannot create ROOT output file: " << outputRootName << endl;
        return;
    }

    int voltage = 0;
    int file_index = 0;
    int selected = 0;
    double amplitude_rpc = 0.0;
    double rise_time_rpc = 0.0;
    double fwhm_rpc = 0.0;
    char csv_path[512] = {0};
    char png_path[512] = {0};

    TTree* tree = new TTree("selection_tree", "RPC waveform shape selection result");
    tree->Branch("voltage", &voltage, "voltage/I");
    tree->Branch("file_index", &file_index, "file_index/I");
    tree->Branch("amplitude_RPC", &amplitude_rpc, "amplitude_RPC/D");
    tree->Branch("rise_time_RPC", &rise_time_rpc, "rise_time_RPC/D");
    tree->Branch("fwhm_RPC", &fwhm_rpc, "fwhm_RPC/D");
    tree->Branch("selected", &selected, "selected/I");
    tree->Branch("csv_path", csv_path, "csv_path/C");
    tree->Branch("png_path", png_path, "png_path/C");

    int total_files = 0;
    int selected_files = 0;

    for (int voltage_value : voltages) {
        voltage = voltage_value;
        const string folderName = to_string(voltage_value);
        const string folderPath = csvFolder + "/" + folderName;
        const int csvCount = countCSVFiles(folderPath);
        const int filesToProcess =
            (maxFilesPerVoltage > 0) ? std::min(csvCount, maxFilesPerVoltage) : csvCount;

        fs::create_directories(fs::path(pngBaseDir) / folderName);

        output->cd();
        TDirectory* waveformDir = output->mkdir(Form("waveforms_%dV", voltage_value));

        cout << "Processing " << voltage_value << "V: " << filesToProcess
             << "/" << csvCount << " CSV files" << endl;

        for (int i = 1; i <= filesToProcess; ++i) {
            file_index = i;
            const string csvFile = folderPath + Form("/%06d.csv", i);
            const string csvStem = Form("%06d", i);
            const string pngFile = (fs::path(pngBaseDir) / folderName / (csvStem + ".png")).string();

            vector<double> time;
            vector<double> rpc_signal;
            vector<double> pmt_signal;

            bool readOk = readCSVData(csvFile, time, rpc_signal, pmt_signal);
            if (!readOk) {
                cerr << "Failed to read " << csvFile << endl;
                continue;
            }

            baselineCorrection(rpc_signal, time);
            baselineCorrection(pmt_signal, time);

            RpcShapeParams params = calculateRpcShapeParams(time, rpc_signal);
            amplitude_rpc = params.amplitude;
            rise_time_rpc = params.rise_time;
            fwhm_rpc = params.fwhm;
            selected = params.valid &&
                       inRange(amplitude_rpc, amp_min, amp_max) &&
                       inRange(rise_time_rpc, rise_min, rise_max) &&
                       inRange(fwhm_rpc, fwhm_min, fwhm_max);

            snprintf(csv_path, sizeof(csv_path), "%s", csvFile.c_str());
            snprintf(png_path, sizeof(png_path), "%s", selected ? pngFile.c_str() : "");
            tree->Fill();
            ++total_files;

            if (!selected) {
                continue;
            }

            ++selected_files;
            TCanvas* canvas = makeWaveformCanvas(voltage_value, csvStem, time,
                                                 rpc_signal, pmt_signal, params);
            canvas->SaveAs(pngFile.c_str());

            waveformDir->cd();
            canvas->Write();
            delete canvas;
            output->cd();
        }
    }

    output->cd();
    tree->Write();
    output->Close();

    cout << "Selection complete." << endl;
    cout << "Total processed CSV files: " << total_files << endl;
    cout << "Selected waveform files: " << selected_files << endl;
    cout << "ROOT output: " << outputRootName << endl;
    cout << "PNG output directory: " << pngBaseDir << endl;
}
