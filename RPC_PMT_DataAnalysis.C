// RPC_PMT_DataAnalysis.C
// 在Root中运行：root -l RPC_PMT_DataAnalysis.C
// 直接处理csv文件，计算信号参数并保存到Root文件中
// 注意： csv文件 三列分别为 [time, rpc, pmt]   调整 readCSVData

// 2026.03.30 by ZhaoYiding 修改
// 优化版本：减少运行时间 2026.04.01 by MIMO
// 修改：增加处理文件数统计，优化进度显示
// 修改：增加TA修正和TQ修正 2026.04.02 by MIMO
// MIMO优化版本v2：减少内存消耗、修复内存泄漏、提升运算速度 2026.04.03

// 2026.04.08 修复画二维图时拟合范围过大导致的拟合失败问题，优化拟合范围设置; 自动获取子文件夹名称并排序
// 2026.04.08 增加：保存TA修正后的各项参数到Root文件中，便于后续分析
// 画图（3*2） 第一行：cfd未修正 TAC TQC；第二行：threshold未修正 TA修正 TQ修正
// 绘制过阈定时、恒比定时的电压依赖曲线，比较修正前后效果（6条曲线画在一幅图上），设置合适的坐标轴范围，添加图例区分不同曲线

// 2026.04.24 时间差采用 PMT-RPC (不取绝对值)；调整TA修正拟合效果

// 2026.05.07 时间差采用 PMT-RPC (不取绝对值)；调整TA修正拟合效果;TAC修正时拟合函数分母没有x
// TAC修正时 对二维直方图在X方向投影，拟合投影得到修正参数，提升拟合稳定性和准确性
// 利用恒比定时计算PMT定时时间，(默认)不采用高斯拟合（循环找点）会节约大量计算时间  //  SignalParams pmt_params = calculateSignalParams(pmt_time, pmt_signal, PMT_AMPLIFICATION, false);
// gaussianFit 如果数据的mean <2 （对应上升时间或半高宽度）, 则sigma初始值设为0.2，否则设为RMS，提升拟合稳定性
/* if (max_bin < 2)
    {
        initial_sigma = 0.2;
    }else if (max_bin >= 1e4 && max_bin < 1e7)
    {
        initial_sigma = 1e5;
    } */

// 创建直方图并进行高斯拟合 处： 幅度bin宽1mV， 电荷bin宽1e3
// 2026.05.09 优化SignalParams calculateSignalParams函数：拟合PMT信号时，每各4个点采集数据，优化初始拟合范围，循环拟合4次提升拟合效果；增加高斯拟合选项参数，默认不使用高斯拟合（循环找点），节约计算时间

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>
#include <map>
#include <limits>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TColor.h"
#include "TProfile.h"
#include "TMath.h"
#include "TDirectory.h"
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include "TProfile.h"
#include <filesystem>
namespace fs = std::filesystem;

using namespace std;

// 全局常量
const double ELECTRON_CHARGE = 1.602176634e-19;
const double IMPEDANCE = 50.0;
const double RPC_AMPLIFICATION = 10.0;
const double PMT_AMPLIFICATION = 1.0;
const double BASELINE_WINDOW = 10.0;
const double PMT_CFD_FRACTION = 0.5;
const double RPC_CFD_FRACTION = 0.2;
const double RPC_THRESHOLD = 10;

const bool FIT_RPC = false; // 是否对RPC信号进行高斯拟合
const bool FIT_PMT = true; // 是否对PMT信号进行高斯拟合

const string csvFolder = "/ustcfs/STCFUser/yzhao/LowPressure/ExpDatas/1090um/5kpa";
const char* outputFileName = "1090um_5kPa";
    

// // vector<int> voltages 
// csvCount = 200;  // 预设文件数量，实际使用时可调用countCSVFiles(FileFolder)获取

// 信号参数结构体
struct SignalParams {
    double amplitude;
    int amplitudex_idx;
    double rise_time;
    double fwhm;
    double charge;
    double timing_50;
    double CFD_time;
    double threshold_time;
};

// 时间差结构体
struct TimingDiff {
    double constant_fraction;
    double threshold;
};

struct FitParameters {
    double mean_amp, sigma_amp;
    double mean_rt, sigma_rt;
    double mean_fwhm, sigma_fwhm;
    double mean_charge, sigma_charge;
    double mean_timing_CFD, sigma_timing_CFD;
    double mean_timing_thresh, sigma_timing_thresh;
};

// 线性插值计算精确时间 - 优化：使用const引用
double interpolateTime(const vector<double>& time, const vector<double>& signal,
                      double threshold, int idx1, int idx2) {
    if (idx1 < 0 || idx2 >= static_cast<int>(signal.size())) return 0.0;

    const double val1 = abs(signal[idx1]);
    const double val2 = abs(signal[idx2]);
    const double t1 = time[idx1];
    const double t2 = time[idx2];

    if (abs(val2 - val1) < 1e-15) return (t1 + t2) * 0.5;

    return t1 + (threshold - val1) * (t2 - t1) / (val2 - val1);
    // return t2;
}

// 函数：计算文件夹中CSV文件的数目
int countCSVFiles(const std::string& folderPath) {
    int count = 0;
    try {
        for (const auto& entry : fs::directory_iterator(folderPath)) {
            if (entry.is_regular_file() && entry.path().extension() == ".csv") {
                count++;
            }
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "错误：无法访问文件夹 - " << e.what() << std::endl;
        return -1;
    }
    return count;
}

// 函数：获取电压列表
vector<int> GetVoltages(const TString& dirPath) {
    vector<int> voltages;

    // 打开目录
    TSystemDirectory dir(dirPath, dirPath);
    TList* files = dir.GetListOfFiles();
    if (!files) {
        cerr << "Cannot open directory: " << dirPath << endl;
        return voltages;
    }

    TIter next(files);
    TSystemFile* file;
    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();

        // 跳过 "." 和 ".."
        if (fname == "." || fname == "..") continue;

        // 只保留子文件夹，且名称为纯数字
        if (file->IsDirectory() && fname.IsDigit()) {
            voltages.push_back(fname.Atoi());
        }
    }

    // 从小到大排序
    sort(voltages.begin(), voltages.end());

    delete files;
    return voltages;
}

// 优化后的文件读取函数 - 使用strtod替代stod提升性能
bool readCSVData(const string& filename, vector<double>& rpc_time, vector<double>& pmt_time,
                         vector<double>& rpc_signal, vector<double>& pmt_signal) {
    ifstream file(filename, ios::binary);
    if (!file.is_open()) {
        return false;
    }

    // 预留内存，减少动态扩容
    rpc_time.reserve(10000);
    pmt_time.reserve(10000);
    rpc_signal.reserve(10000);
    pmt_signal.reserve(10000);

    string line;
    // getline(file, line); // 跳过标题行

    // 输出前10行数据进行检查
    int line_count = 0;

    // time, rpc, pmt
    // 注意：根据实际数据格式调整读取顺序和符号
    while (getline(file, line)) {
        if (line.empty()) continue;

        // 使用strtod替代stod，性能更优且无需异常处理
        const char* ptr = line.c_str();
        char* end_ptr;

        double t = strtod(ptr, &end_ptr) * 1e9;
        if (end_ptr == ptr) continue;

        ptr = end_ptr + 1;  // 跳过逗号
        double rpc = -strtod(ptr, &end_ptr) * 1e3; //负信号变为正信号
        if (end_ptr == ptr) continue;

        ptr = end_ptr + 1;  // 跳过逗号
        double pmt = -strtod(ptr, &end_ptr) * 1e3; //负信号变为正信号
        if (end_ptr == ptr) continue;

        // 全部换成正信号
        rpc_time.push_back(t);
        pmt_time.push_back(t);
        rpc_signal.push_back(rpc);
        pmt_signal.push_back(pmt);

        // // 输出前10行数据进行检查
        // if (line_count < 10) {
        //     cout << "Line " << line_count + 1 << ": time = " << t << " ns, PMT = " << pmt << " mV, RPC = " << rpc << " mV" << endl;
        // }
        // line_count++;

    }

    file.close();
    if (rpc_time.size() < 3 || rpc_time.size() != pmt_time.size() ||
        rpc_time.size() != rpc_signal.size() || rpc_time.size() != pmt_signal.size()) return false;
    for (size_t i = 1; i < rpc_time.size(); ++i) {
        if (!std::isfinite(rpc_time[i]) || rpc_time[i] <= rpc_time[i - 1]) return false;
    }
    return true;
}

// 优化的基线校正 - 预存size避免重复调用
void baselineCorrection(vector<double>& signal, const vector<double>& time, const double baseline_window = BASELINE_WINDOW) {
    if (signal.size() < 3 || signal.size() != time.size()) return;

    const double step = time[2] - time[1];
    const double baseline_window_end = time[0] + baseline_window;
    const int baseline_window_end_idx = static_cast<int>((baseline_window_end - time[0]) / step);
    const int sig_size = static_cast<int>(signal.size());

    double baseline_sum = 0.0;
    int baseline_count = 0;

    for (int i = 0; i < baseline_window_end_idx && i < sig_size; ++i) {
        baseline_sum += signal[i];
        baseline_count++;
    }

    if (baseline_count > 0) {
        const double baseline_mean = baseline_sum / baseline_count;
        for (size_t i = 0; i < signal.size(); ++i) {
            signal[i] -= baseline_mean;
        }
    }
}

// 优化的信号参数计算 - 修复内存泄漏，优化循环
SignalParams calculateSignalParams( const vector<double>& time,
                                    const vector<double>& signal,
                                    double amplification_factor, 
                                    bool Gaus_Fit = false) {
    SignalParams params = {0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    if (signal.size() < 3 || signal.size() != time.size()) return params;

    // const double time_min = -8.0;
    const double time_min = max(-14.0, time[0]);
    const double time_max = min(14.0, time.back());
    const double step = time[2] - time[1];
    int time_min_idx = static_cast<int>((time_min - time[0]) / step);
    int time_max_idx = static_cast<int>((time_max - time[0]) / step);

    time_min_idx = max(0, time_min_idx);
    time_max_idx = min(static_cast<int>(signal.size()) - 1, time_max_idx);

    auto amp_max_it = max_element(signal.begin() + time_min_idx, signal.begin() + time_max_idx + 1);
    auto amp_min_it = min_element(signal.begin() + time_min_idx, signal.begin() + time_max_idx + 1);

    const double max_amp_val = *amp_max_it;
    const double min_amp_val = *amp_min_it;
    const int max_amp_idx = amp_max_it - signal.begin();
    const int min_amp_idx = amp_min_it - signal.begin();

    params.amplitude = max(abs(max_amp_val), abs(min_amp_val));
    params.amplitudex_idx = abs(max_amp_val) > abs(min_amp_val) ? max_amp_idx : min_amp_idx;

    // 对PMT信号进行高斯拟合，获取更准确的幅值和时间
    double A_f_pmt = 0;
    double mean_f_pmt = 0;
    double sigma_f_pmt = 0;
    bool gaussian_fit_valid = false;

    if ( Gaus_Fit ){

        // zyd
        // const double Range_min = - 8; // 起始时间
        // const double Range_max = 8; // 结束时间
        const double Range_min = time[params.amplitudex_idx]-8; // 起始时间
        const double Range_max = time[params.amplitudex_idx]+8; // 结束时间


        // 1. 过滤数据
        std::vector<double> filtered_time, filtered_signal;
        for (size_t i = 0; i < time.size(); ++i) {
            if (time[i] >= Range_min && time[i] <= Range_max) {
                filtered_time.push_back(time[i]);
                filtered_signal.push_back(signal[i]);
            }
        }

        const double t1 = time[params.amplitudex_idx] - 4; // 起始时间
        const double t2 = time[params.amplitudex_idx] + 1; // 结束时间

        // 2. 创建 TGraphErrors（假设误差为 0，如有实际误差可添加）
        // TGraphErrors *g_pmt = new TGraphErrors(time.size());
        // for (size_t i = 0; i < time.size(); ++i) {
        //     g_pmt->SetPoint(i, filtered_time[i], filtered_signal[i]);
        //     g_pmt->SetPointError(i, 0, 0); // 如果没有误差，设为 0
        // }

        // 2. 创建 TGraph（无误差，比TGraphErrors更轻量）
        TGraph *g_pmt = new TGraph(filtered_time.size());
        for (size_t i = 0; i < filtered_time.size(); i=i+1) {
            g_pmt->SetPoint(i, filtered_time[i], filtered_signal[i]);
        }

        TF1 *f_pmt = new TF1("f_pmt", "gaus",t1,t2);
        f_pmt->SetParameters(signal[params.amplitudex_idx], time[params.amplitudex_idx], 1.7); // 初始参数：幅值、均值、标准差
        int fit_status = g_pmt->Fit("f_pmt", "R Q");
        // Get the updated fit results
        A_f_pmt = f_pmt->GetParameter(0);
        mean_f_pmt = f_pmt->GetParameter(1);
        sigma_f_pmt = f_pmt->GetParameter(2);
        // cout<<"Initial PMT Fit: A = " << A_f_pmt << ", mean = " << mean_f_pmt << " ns, sigma = " << sigma_f_pmt << " ns" << endl;
        // cout.precision(6);

        // cout<<"A_f_pmt, mean_f_pmt, sigma_f_pmt: " << A_f_pmt << ", " << mean_f_pmt << " ns, " << sigma_f_pmt << " ns" << endl;

        // double PMT_Fit_T1 = 0.1;
        // double PMT_Fit_T2 = 0.7;
        double PMT_Fit_T1 = 0.1;
        double PMT_Fit_T2 = 0.7;
        double P_or_N = +1;

        // Loop fit twice
        // zyd 大量数据多次拟合比较耗时间； 
        // 不循环拟合计算很快，但有些数据拟合效果不好；循环两次拟合可以提升拟合效果，但会增加运行时间
        for (int i = 0; i < 2; i++)
        {
            const double fit_t1 = -sqrt(-log(PMT_Fit_T1))*sqrt(2)*sigma_f_pmt + mean_f_pmt;
            const double fit_t2 = P_or_N * sqrt(-log(PMT_Fit_T2)) * sqrt(2) * sigma_f_pmt + mean_f_pmt;
            f_pmt->SetRange(fit_t1, fit_t2);
            f_pmt->SetParameters(A_f_pmt, mean_f_pmt, sigma_f_pmt);

            // cout<<"Fit iteration " << i+1 << ": fit_t1 = " << fit_t1 << " ns, fit_t2 = " << fit_t2 << " ns" << endl;    

            // refit
            fit_status = g_pmt->Fit(f_pmt,"R Q");
            // Get the updated fit results
            A_f_pmt = f_pmt->GetParameter(0);
            mean_f_pmt = f_pmt->GetParameter(1);
            sigma_f_pmt = f_pmt->GetParameter(2);

            // cout<<"A_f_pmt, mean_f_pmt, sigma_f_pmt: " << A_f_pmt << ", " << mean_f_pmt << " ns, " << sigma_f_pmt << " ns" << endl;
        }
        gaussian_fit_valid = fit_status == 0 && std::isfinite(A_f_pmt) &&
                             std::isfinite(mean_f_pmt) && std::isfinite(sigma_f_pmt) &&
                             sigma_f_pmt > 0.0;
        // cout<<"\t"<<endl;
    }
  
    const double threshold_10 = 0.1 * params.amplitude;
    const double threshold_90 = 0.9 * params.amplitude;
    const double threshold_50 = 0.5 * params.amplitude;
    const double threshold_01 = 0.01 * params.amplitude;
    const double threshold_CFD = RPC_CFD_FRACTION * params.amplitude;
    const double threshold_A = RPC_THRESHOLD;

    int idx_10_start = -1, idx_90_start = -1, idx_50_left = -1;
    int idx_01_left = -1, idx_CFD_start = -1, idx_threshold_start = -1, idx_threshold_prev = -1;
    int idx_50_right = -1, idx_01_right = -1;

    const int start_idx = params.amplitudex_idx;
    const int sig_size = static_cast<int>(signal.size());

    // 优化：预计算所有abs值避免重复计算
    double integral_left = 0.0;
    double integral_right = 0.0;
    for (int i = start_idx; i > 0; --i) {
        const double val = abs(signal[i]);  // 只计算一次abs

        if (idx_90_start == -1 && val <= threshold_90) idx_90_start = i;
        if (idx_50_left == -1 && val <= threshold_50) idx_50_left = i;
        if (idx_CFD_start == -1 && val <= threshold_CFD) idx_CFD_start = i;
        if (idx_10_start == -1 && val <= threshold_10) idx_10_start = i;
        if (idx_01_left == -1 && val <= threshold_01) idx_01_left = i;

        // double dt = (time[i] - time[i-1]) * 1e-9;
        const double dt = (time[i] - time[i - 1]) * 1e-9;
        const double current = 0.5 * (signal[i] + signal[i - 1]) * 1e-3 / IMPEDANCE;
        integral_left += current * dt;

        if (idx_01_left != -1 && idx_90_start != -1 && idx_10_start != -1 &&
            idx_50_left != -1 && idx_CFD_start != -1) break;
    }

    int threshold_search_start_idx = lower_bound(time.begin(), time.end(), -9.0) - time.begin();
    threshold_search_start_idx = max(1, min(threshold_search_start_idx, sig_size - 1));
    for (int i = threshold_search_start_idx; i < sig_size; ++i) {
        if (abs(signal[i]) >= threshold_A) {
            idx_threshold_start = i;
            idx_threshold_prev = i - 1;
            break;
        }
    }

    for (int i = start_idx; i + 1 < sig_size; ++i) {
        const double val = abs(signal[i]);  // 只计算一次abs

        if (idx_50_right == -1 && val <= threshold_50) idx_50_right = i;
        if ((idx_01_right == -1 && val <= threshold_01) || (i > sig_size - 2))
            idx_01_right = i;

        // const double dt = (time[i+1] - time[i]) * 1e-9;
        const double dt = (time[i + 1] - time[i]) * 1e-9;
        const double current = 0.5 * (signal[i] + signal[i + 1]) * 1e-3 / IMPEDANCE;
        integral_right += current * dt;

        if (idx_01_right != -1 && idx_50_right != -1) break;
    }

    if (idx_10_start > 0 && idx_90_start > 0) {
        const double t10 = interpolateTime(time, signal, threshold_10, idx_10_start - 1, idx_10_start);
        const double t90 = interpolateTime(time, signal, threshold_90, idx_90_start - 1, idx_90_start);
        params.rise_time = abs(t90 - t10);
    }

    if (idx_50_left > 0 && idx_50_right > 0) {
        const double t1 = interpolateTime(time, signal, threshold_50, idx_50_left - 1, idx_50_left);
        const double t2 = interpolateTime(time, signal, threshold_50, idx_50_right + 1, idx_50_right);
        params.fwhm = abs(t2 - t1);
    }

    // //zyd
    // cout<<"amplitude: " << params.amplitude << " mV, rise_time: " << params.rise_time << " ns, FWHM: " << params.fwhm << endl;
    // cout<<"time_10  "<< time[idx_10_start] << " ns, time_50_left: " << time[idx_50_left] << " ns, time_50_right: " << time[idx_50_right] << " ns" << endl;
    // cout<<"time_90  "<< time[idx_90_start] << " ns, time_CFD: " << time[idx_CFD_start] << " ns, time_threshold: " << time[idx_threshold_start] << " ns" << endl<<endl;

    // cout<<"amplitude: "<<signal[params.amplitudex_idx]<<endl;
    // cout<<"am_10 " <<signal[idx_10_start]<< "  Ratio = "<<signal[idx_10_start]/params.amplitude<<endl;
    // cout<<"am_50_left " <<signal[idx_50_left]<<"  Ratio = "<<signal[idx_50_left]/params.amplitude<< endl;
    // cout<<"am_50_right " <<signal[idx_50_right]<<"  Ratio = "<<signal[idx_50_right]/params.amplitude<< endl;
    // cout<<"am_90 " <<signal[idx_90_start]<<"  Ratio = "<<signal[idx_90_start]/params.amplitude<< endl<<endl;;

    if (idx_01_left != -1 && idx_01_right != -1 && idx_01_right > idx_01_left) {
        double integral = 0.0;
        // for (int i = idx_01_left; i < idx_01_right; ++i) {
        //     const double dt = (time[i+1] - time[i]) * 1e-9;
        //     const double current = (signal[i] * 1e-3) / IMPEDANCE;
        //     integral += current * dt;
        // }
        integral = integral_left + integral_right;
        params.charge = abs(integral) / ELECTRON_CHARGE / amplification_factor;
    }

    params.timing_50 = interpolateTime(time, signal, threshold_50, idx_50_left - 1, idx_50_left);
    params.CFD_time = interpolateTime(time, signal, threshold_CFD, idx_CFD_start - 1, idx_CFD_start);
    params.threshold_time = RPC_THRESHOLD < params.amplitude ?
                           interpolateTime(time, signal, threshold_A, idx_threshold_prev, idx_threshold_start) : std::numeric_limits<double>::quiet_NaN();

    if (Gaus_Fit && gaussian_fit_valid){

        const double pmt_time_50 = -sqrt(-log(PMT_CFD_FRACTION))*sqrt(2)*sigma_f_pmt + mean_f_pmt;
        const double pmt_time_CFD = -sqrt(-log(RPC_CFD_FRACTION))*sqrt(2)*sigma_f_pmt + mean_f_pmt;
        const double pmt_time_thresh = -sqrt(-log(RPC_THRESHOLD/params.amplitude))*sqrt(2)*sigma_f_pmt + mean_f_pmt;

        params.timing_50 = pmt_time_50;
        params.CFD_time = pmt_time_CFD;
        params.threshold_time = pmt_time_thresh;
    }

    return params;
}

// 时间差计算
TimingDiff calculateTimingDiff(const SignalParams& rpc_params,
                              const SignalParams& pmt_params) {
    TimingDiff diff = {0.0, 0.0};
    // diff.constant_fraction = abs(rpc_params.CFD_time - pmt_params.timing_50);
    // diff.threshold = abs(rpc_params.threshold_time - pmt_params.timing_50);
    diff.constant_fraction = -(rpc_params.CFD_time - pmt_params.timing_50);
    diff.threshold = -(rpc_params.threshold_time - pmt_params.timing_50);
    return diff;
}

// 高斯拟合
void gaussianFit(TH1F* hist, double& mean, double& sigma) {
    if (!hist || hist->GetEntries() == 0) return;

    const double max_content = hist->GetMaximum();
    const double max_bin = hist->GetBinCenter(hist->GetMaximumBin());
    double initial_sigma = hist->GetRMS();

    if (max_bin < 2)
    {
        initial_sigma = 0.2;
    }else if (max_bin >= 1e4 && max_bin < 1e7)
    {
        initial_sigma = 1e5;
    }
    
    

    TF1* gauss = new TF1("gauss", "gaus", max_bin - 1*initial_sigma, max_bin + 1*initial_sigma);
    gauss->SetParameters(max_content, max_bin, initial_sigma);
    hist->Fit(gauss, "RQ");
    gStyle->SetOptFit(111);

    mean = gauss->GetParameter(1);
    sigma = gauss->GetParameter(2);

    for (int i = 0; i < 2; ++i) {
        gauss->SetRange(mean - 2*sigma, mean + 2*sigma);
        hist->Fit(gauss, "RQ");
        mean = gauss->GetParameter(1);
        sigma = gauss->GetParameter(2);
    }

    delete gauss;
}

// ========== 新增：迭代高斯拟合函数 ==========
void iterativeGaussianFit(  TH1F* hist, int nIter,
                            double& mean, double& sigma,
                            double& meanErr, double& sigmaErr,
                            double sigmaScale = 3.0) {
    mean = 0; sigma = 0; meanErr = 0; sigmaErr = 0;
    if (!hist || hist->GetEntries() == 0) return;

    const double max_content = hist->GetMaximum();
    const double max_bin = hist->GetBinCenter(hist->GetMaximumBin());
    const double initial_sigma = hist->GetRMS();

    // // 创建独立的TF1对象，使用空名称让ROOT自动命名
    // TF1* gfit = new TF1("", "gaus",
    //                     hist->GetXaxis()->GetXmin(),
    //                     hist->GetXaxis()->GetXmax());
    // hist->Fit(gfit, "RQ");
    // gStyle->SetOptFit(111);

    TF1* gfit = new TF1("gauss", "gaus", max_bin - 2*initial_sigma, max_bin + 2*initial_sigma);
    gfit->SetParameters(max_content, max_bin, initial_sigma);
    hist->Fit(gfit, "RQ");
    gStyle->SetOptFit(111);

    mean = gfit->GetParameter(1);
    sigma = gfit->GetParameter(2);
    meanErr = gfit->GetParError(1);
    sigmaErr = gfit->GetParError(2);

    for (int i = 0; i < nIter; ++i) {
        const double fitRange = 1 * sigma;
        const double lo = mean - sigmaScale * fitRange;
        const double hi = mean + sigmaScale * fitRange;
        if (lo < hi) {
            gfit->SetRange(lo, hi);
            hist->Fit(gfit, "RQ");
            mean = gfit->GetParameter(1);
            sigma = gfit->GetParameter(2);
            meanErr = gfit->GetParError(1);
            sigmaErr = gfit->GetParError(2);
        }
    }
    // 不删除gfit，因为它被hist内部引用
}

// 串行处理电压的函数（修复版）
void processVoltageSerial(int voltage, const string& folder_name,
                         vector<double>& amplitudes,
                         vector<double>& rise_times,
                         vector<double>& fwhms,
                         vector<double>& charges,
                         vector<double>& timing_diffs_const,
                         vector<double>& timing_diffs_thresh,
                         TTree* tree) {

    // string csvFolder = "/ustcfs/STCFUser/yzhao/photo_RPC/Data_2026/e10_1atm_STDgas";
    
    string FileFolder = csvFolder + Form("/%s/", folder_name.c_str());
    int csvCount = countCSVFiles(FileFolder);
    // int csvCount = 2;  // 预设文件数量，实际使用时可调用countCSVFiles(FileFolder)获取

    // 【修复】分支变量声明在循环外部，确保生命周期覆盖所有Fill调用
    double amp_rpc, amp_pmt, rt_rpc, rt_pmt, fwhm_rpc, fwhm_pmt;
    double charge_rpc, charge_pmt, timing_CFD, timing_thresh;
    bool threshold_valid = false;

    // 【修复】使用Branch而非SetBranchAddress（只设置一次指针）
    tree->Branch("amplitude_RPC", &amp_rpc, "amplitude_RPC/D");
    tree->Branch("amplitude_PMT", &amp_pmt, "amplitude_PMT/D");
    tree->Branch("rise_time_RPC", &rt_rpc, "rise_time_RPC/D");
    tree->Branch("rise_time_PMT", &rt_pmt, "rise_time_PMT/D");
    tree->Branch("fwhm_RPC", &fwhm_rpc, "fwhm_RPC/D");
    tree->Branch("fwhm_PMT", &fwhm_pmt, "fwhm_PMT/D");
    tree->Branch("charge_RPC", &charge_rpc, "charge_RPC/D");
    tree->Branch("charge_PMT", &charge_pmt, "charge_PMT/D");
    tree->Branch("timing_CFD", &timing_CFD, "timing_CFD/D");
    tree->Branch("timing_threshold", &timing_thresh, "timing_threshold/D");
    tree->Branch("threshold_valid", &threshold_valid, "threshold_valid/O");

    // 预留容器空间减少动态扩容
    amplitudes.reserve(csvCount);
    rise_times.reserve(csvCount);
    fwhms.reserve(csvCount);
    charges.reserve(csvCount);
    timing_diffs_const.reserve(csvCount);
    timing_diffs_thresh.reserve(csvCount);

    // for (int file_num = 1; file_num <= 300; ++file_num) {
    for (int file_num = 1; file_num <= csvCount; ++file_num) {

        if (file_num % 50 == 0) {
            cout << "folder_name: " << folder_name << "V, Processing file: "
                 << file_num << "/" << csvCount
                 << "  ****progress**** " << file_num * 100 / csvCount << "%" << endl;
        }

        string filename = FileFolder + Form("/%06d.csv", file_num);

        vector<double> rpc_time, pmt_time, rpc_signal, pmt_signal;
        if (!readCSVData(filename, rpc_time, pmt_time, rpc_signal, pmt_signal)) {
            if (file_num == 1) {
                cout << "文件夹 " << folder_name << " 中没有找到数据文件" << endl;
            }
            break;
        }

        baselineCorrection(rpc_signal, rpc_time);
        baselineCorrection(pmt_signal, pmt_time);

        // zyd
        SignalParams rpc_params = calculateSignalParams(rpc_time, rpc_signal, RPC_AMPLIFICATION, FIT_RPC);
        // SignalParams pmt_params = calculateSignalParams(pmt_time, pmt_signal, PMT_AMPLIFICATION, true);
        SignalParams pmt_params = calculateSignalParams(pmt_time, pmt_signal, PMT_AMPLIFICATION, FIT_PMT);

        TimingDiff timing_diff = calculateTimingDiff(rpc_params, pmt_params);

        amp_rpc = rpc_params.amplitude;
        amp_pmt = pmt_params.amplitude;
        rt_rpc = rpc_params.rise_time;
        rt_pmt = pmt_params.rise_time;
        fwhm_rpc = rpc_params.fwhm;
        fwhm_pmt = pmt_params.fwhm;
        charge_rpc = rpc_params.charge;
        charge_pmt = pmt_params.charge;
        timing_CFD = timing_diff.constant_fraction;
        timing_thresh = timing_diff.threshold;
        threshold_valid = std::isfinite(timing_thresh);

        tree->Fill();

        amplitudes.push_back(amp_rpc);
        rise_times.push_back(rt_rpc);
        fwhms.push_back(fwhm_rpc);
        charges.push_back(charge_rpc);
        timing_diffs_const.push_back(timing_CFD);
        timing_diffs_thresh.push_back(timing_thresh);
    }

    // 【修复】在所有文件处理完毕后才写入TTree
    tree->Write();
}



// ========== 修改后的 performTA_TQ_Correction ==========
// 对每个电压：
//   1. 创建独立目录 Cali_XXXV
//   2. 图1(2x2)：幅度分布、电荷分布、CFD时间差分布、阈值时间差分布，各自高斯拟合
//   3. 图2(2x2)：四个二维分布 + 拟合曲线
//   4. 图3(2x2)：修正前后时间差分布对比
//   5. 汇总：分辨率 vs 电压图（一级目录 + PNG）
// 优化：修复内存泄漏，预存size
void performTA_TQ_Correction(const vector<int>& voltages,
                             const map<int, FitParameters>& fitParamsMap,
                             const map<int, vector<double>>& amplitudes,
                             const map<int, vector<double>>& charges,
                             const map<int, vector<double>>& rise_times,
                             const map<int, vector<double>>& fwhms,
                             const map<int, vector<double>>& timing_diffs_const,
                             const map<int, vector<double>>& timing_diffs_thresh,
                             TFile* outputFile,
                             map<int, double>& cfd_sigma_before,
                             map<int, double>& cfd_sigmaErr_before,
                             map<int, double>& cfd_sigma_afterTA,
                             map<int, double>& cfd_sigmaErr_afterTA,
                             map<int, double>& cfd_sigma_afterTQ,
                             map<int, double>& cfd_sigmaErr_afterTQ,
                             map<int, double>& thresh_sigma_before,
                             map<int, double>& thresh_sigmaErr_before,
                             map<int, double>& thresh_sigma_afterTA,
                             map<int, double>& thresh_sigmaErr_afterTA,
                             map<int, double>& thresh_sigma_afterTQ,
                             map<int, double>& thresh_sigmaErr_afterTQ) {

    gStyle->SetOptFit(111);
    for (int voltage : voltages) {
        gStyle->SetOptFit(111);

        cout << "\n========== TA/TQ 修正: " << voltage << "V ==========" << endl;

        const auto& amps0 = amplitudes.at(voltage);
        const auto& chgs0 = charges.at(voltage);
        const auto& tCFD0 = timing_diffs_const.at(voltage);
        const auto& tThresh0 = timing_diffs_thresh.at(voltage);
        const auto& rts0 = rise_times.at(voltage);
        const auto& fwhm0 = fwhms.at(voltage); // 使用拟合的FWHM参数进行筛选

        vector<double> amps, chgs, tCFD, tThresh, rts;
        // 预留空间减少动态扩容
        amps.reserve(amps0.size());
        chgs.reserve(chgs0.size());
        tCFD.reserve(tCFD0.size());
        tThresh.reserve(tThresh0.size());
        rts.reserve(rts0.size());

        const double RT_min = fitParamsMap.at(voltage).mean_rt - 5 * fitParamsMap.at(voltage).sigma_rt;
        const double RT_max = fitParamsMap.at(voltage).mean_rt + 5 * fitParamsMap.at(voltage).sigma_rt;

        const double A_min = fitParamsMap.at(voltage).mean_amp - 5 * fitParamsMap.at(voltage).sigma_amp;
        const double A_max = fitParamsMap.at(voltage).mean_amp + 5 * fitParamsMap.at(voltage).sigma_amp;

        const double fwhm_min = fitParamsMap.at(voltage).mean_fwhm - 5 * fitParamsMap.at(voltage).sigma_fwhm;
        const double fwhm_max = fitParamsMap.at(voltage).mean_fwhm + 5 * fitParamsMap.at(voltage).sigma_fwhm;

        const int data_size = static_cast<int>(amps0.size());

        for (int i = 0; i < data_size; ++i) {
            // 筛选：使用 amps0[i] 和 rts0[i] 进行条件检查
            const bool A_cut = (amps0[i] > A_min && amps0[i] < A_max);  // 修正：使用 amps0[i]
            const bool RT_cut = (rts0[i] > RT_min && rts0[i] < RT_max);  // 修正：使用 rts0[i]
            const bool FWHM_cut = (fwhm0[i] > fwhm_min && fwhm0[i] < fwhm_max);  // 修正：使用 fwhm0[i]

            // if (A_cut && RT_cut) {
            // if (RT_cut && FWHM_cut) {  // 优化：去掉幅度筛选，保留时间相关筛选
             if (RT_cut ) {  // 优化：去掉幅度筛选，保留时间相关筛选
                amps.push_back(amps0[i]);
                chgs.push_back(chgs0[i]);
                tCFD.push_back(tCFD0[i]);
                tThresh.push_back(tThresh0[i]);
                rts.push_back(rts0[i]);
            }
        }


        // ===== 图2：四个2D分布 + 拟合曲线 (2x2) ==========

        const int n = static_cast<int>(amps.size());
        if (n < 10) {
            cout << "  数据点不足，跳过" << endl;
            continue;
        }

        // ===== 创建该电压的独立目录 =====
        TDirectory* caliDir = outputFile->mkdir(Form("Cali_%dV", voltage));
        caliDir->cd();

        // ===== 计算统计量 =====
        const double meanA = fitParamsMap.at(voltage).mean_amp;
        const double rmsA  = fitParamsMap.at(voltage).sigma_amp;
        const double meanQ = fitParamsMap.at(voltage).mean_charge;
        const double rmsQ  = fitParamsMap.at(voltage).sigma_charge;
        const double meanTcfd = fitParamsMap.at(voltage).mean_timing_CFD;
        const double rmsTcfd  = fitParamsMap.at(voltage).sigma_timing_CFD;

        // 阈值定时差：筛选有效值
        const double meanT_thresh = fitParamsMap.at(voltage).mean_timing_thresh;
        const double rmsT_thresh  = fitParamsMap.at(voltage).sigma_timing_thresh;


        vector<double> tThreshValid;
        tThreshValid.reserve(n);
        for (int i = 0; i < n; ++i) {
            if (std::isfinite(tThresh[i])) tThreshValid.push_back(tThresh[i]);
        }

        // ===== 1D直方图范围 =====
        const int sigma_num = 10;
        const double A_lo  = max(0.1, meanA - sigma_num*rmsA);
        const double A_hi  = meanA + sigma_num*rmsA;
        const double Q_lo  = max(0.1, meanQ - sigma_num*rmsQ);
        const double Q_hi  = meanQ + sigma_num*rmsQ;
        const double Tcfd_lo = meanTcfd - sigma_num*rmsTcfd;
        const double Tcfd_hi = meanTcfd + sigma_num*rmsTcfd;
        const double Tthresh_lo = meanT_thresh - sigma_num*rmsT_thresh;
        const double Tthresh_hi = meanT_thresh + sigma_num*rmsT_thresh;


        // zyd
        const int nBinsA = max(100, (int)((A_hi - A_lo) / max(rmsA/20.0, 0.5)));
        const int nBinsC = max(100, (int)((Q_hi - Q_lo) / max(rmsQ/20.0, 1e3)));
        const double binW_Tcfd = 0.005; // 0.005
        const double binW_Tthresh = 0.005;
        const int nBinsTcfd = max(200, (int)((Tcfd_hi - Tcfd_lo) / binW_Tcfd));
        const int nBinsTthresh = max(200, (int)((Tthresh_hi - Tthresh_lo) / binW_Tthresh));

        // ===== 创建1D直方图 =====
        TH1F* hAmp = new TH1F(Form("hF1_amp_%dV", voltage),
            Form("Amplitude at %dV;Amplitude (mV);Counts", voltage), nBinsA, A_lo, A_hi);
        TH1F* hChg = new TH1F(Form("hF1_chg_%dV", voltage),
            Form("Charge at %dV;Charge (Q_{e});Counts", voltage), nBinsC, Q_lo, Q_hi);
        TH1F* hTC = new TH1F(Form("hF1_tc_%dV", voltage),
            Form("CFD Timing at %dV;Time Difference (ns);Counts", voltage), nBinsTcfd, Tcfd_lo, Tcfd_hi);
        TH1F* hTT = new TH1F(Form("hF1_tt_%dV", voltage),
            Form("Threshold Timing at %dV;Time Difference (ns);Counts", voltage), nBinsTthresh, Tthresh_lo, Tthresh_hi);

        for (int i = 0; i < n; ++i) {
            hAmp->Fill(amps[i]);
            hChg->Fill(chgs[i]);
            hTC->Fill(tCFD[i]);
            if (std::isfinite(tThresh[i])) hTT->Fill(tThresh[i]);
        }

        // ========== 图1：四个1D分布 (2x2) ==========
        TCanvas* cFig1 = new TCanvas(Form("cF1_%dV", voltage),
            Form("Distributions at %dV", voltage), 1600, 1200);
        cFig1->Divide(2, 2);

        // Pad 1: 幅度分布
        cFig1->cd(1); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        hAmp->Draw();
        {
            double m, s, mErr, sErr;
            iterativeGaussianFit(hAmp, 2, m, s, mErr, sErr);

            // // 重新设置拟合范围为±2σ
            // TF1* f = new TF1("", "gaus", m - 2*s, m + 2*s);
            // hAmp->Fit(f, "RQ");
            // m = f->GetParameter(1); s = f->GetParameter(2);
            // f->SetRange(m - 2*s, m + 2*s);
            // hAmp->Fit(f, "RQ");
            
            // delete f;  // 释放内存
        }

        // Pad 2: 电荷分布
        cFig1->cd(2); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        hChg->Draw();
        {
            double m, s, mErr, sErr;
            iterativeGaussianFit(hChg, 2, m, s, mErr, sErr);
            // TF1* f = new TF1("", "gaus", m - 2*s, m + 2*s);
            // hChg->Fit(f, "RQ");
            // m = f->GetParameter(1); s = f->GetParameter(2);
            // f->SetRange(m - 2*s, m + 2*s);
            // hChg->Fit(f, "RQ");
            // delete f;  // 释放内存
        }

        // Pad 3: CFD时间差分布
        cFig1->cd(3); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        hTC->Draw();
        {
            double m, s, mErr, sErr;
            iterativeGaussianFit(hTC, 2, m, s, mErr, sErr);
            // TF1* f = new TF1("", "gaus", m - 2*s, m + 2*s);
            // hTC->Fit(f, "RQ");
            // m = f->GetParameter(1); s = f->GetParameter(2);
            // f->SetRange(m - 2*s, m + 2*s);
            // hTC->Fit(f, "RQ");
            // delete f;  // 释放内存
        }

        // Pad 4: 阈值时间差分布
        cFig1->cd(4); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        hTT->Draw();
        {
            double m, s, mErr, sErr;
            iterativeGaussianFit(hTT, 2, m, s, mErr, sErr);
            // TF1* f = new TF1("", "gaus", m - 2*s, m + 2*s);
            // hTT->Fit(f, "RQ");
            // m = f->GetParameter(1); s = f->GetParameter(2);
            // f->SetRange(m - 2*s, m + 2*s);
            // hTT->Fit(f, "RQ");
            // delete f;  // 释放内存
        }

        cFig1->Write();
        // cFig1->SaveAs(Form("Cali_%dV_Fig1_Distributions.png", voltage));

        // ========== 图2：四个2D分布 + 拟合曲线 (2x2) ==========

        const int sigma_left = 4;
        const int sigma_right = 4;
        const double h2A_lo = max(6.0, meanA - sigma_left*rmsA);
        const double h2A_hi = meanA + sigma_right*rmsA;
        const double h2Q_lo = max(2e3, meanQ - sigma_left*rmsQ);
        const double h2Q_hi = meanQ + sigma_right*rmsQ;

        // const double h2Tcfd_lo = max(0.1, meanTcfd - sigma_left*rmsTcfd);
        // const double h2Tcfd_hi = meanTcfd + sigma_right*rmsTcfd;
        // const double h2Tthresh_lo = max(0.1, meanT_thresh - sigma_left*rmsT_thresh);
        // const double h2Tthresh_hi = meanT_thresh + sigma_right*rmsT_thresh;

        const double h2Tcfd_lo =  meanTcfd - sigma_left*rmsTcfd;
        const double h2Tcfd_hi = meanTcfd + sigma_right*rmsTcfd;
        const double h2Tthresh_lo =  meanT_thresh - sigma_left*rmsT_thresh;
        const double h2Tthresh_hi = meanT_thresh + sigma_right*rmsT_thresh;

        double bin_num = 100;

        TH2F* h2_CFD_A = new TH2F(Form("hF2_cfdA_%dV", voltage),
            Form("CFD vs Amplitude %dV;Amplitude (mV);Time Diff (ns)", voltage),
            bin_num, h2A_lo, h2A_hi, bin_num, h2Tcfd_lo, h2Tcfd_hi);
        TH2F* h2_CFD_Q = new TH2F(Form("hF2_cfdQ_%dV", voltage),
            Form("CFD vs Charge %dV;Charge (Q_{e});Time Diff (ns)", voltage),
            bin_num, h2Q_lo, h2Q_hi, bin_num, h2Tcfd_lo, h2Tcfd_hi);
        TH2F* h2_Thresh_A = new TH2F(Form("hF2_thA_%dV", voltage),
            Form("Threshold vs Amplitude %dV;Amplitude (mV);Time Diff (ns)", voltage),
            bin_num, h2A_lo, h2A_hi, bin_num, h2Tthresh_lo, h2Tthresh_hi);
        TH2F* h2_Thresh_Q = new TH2F(Form("hF2_thQ_%dV", voltage),
            Form("Threshold vs Charge %dV;Charge (Q_{e});Time Diff (ns)", voltage),
            bin_num, h2Q_lo, h2Q_hi, bin_num, h2Tthresh_lo, h2Tthresh_hi);

        for (int i = 0; i < n; ++i) {
            h2_CFD_A->Fill(amps[i], tCFD[i]);
            h2_CFD_Q->Fill(chgs[i], tCFD[i]);
            if (std::isfinite(tThresh[i])) {
                h2_Thresh_A->Fill(amps[i], tThresh[i]);
                h2_Thresh_Q->Fill(chgs[i], tThresh[i]);
            }


        }


        // 二维直方图向X方向投影得到TProfile
        TProfile *hCFD_A_ProfileX = h2_CFD_A->ProfileX();
        TProfile *hCFD_Q_ProfileX = h2_CFD_Q->ProfileX();
        TProfile *hThresh_A_ProfileX = h2_Thresh_A->ProfileX();
        TProfile *hThresh_Q_ProfileX = h2_Thresh_Q->ProfileX();

        

        // 构建 TGraph 用于拟合
        TGraph* gr_cfdA = new TGraph(n, amps.data(), tCFD.data());
        TGraph* gr_cfdQ = new TGraph(n, chgs.data(), tCFD.data());
        const int tThreshValidSize = static_cast<int>(tThreshValid.size());
        TGraph* gr_thA  = new TGraph(tThreshValidSize);
        TGraph* gr_thQ  = new TGraph(tThreshValidSize);
        {
            int idx = 0;
            for (int i = 0; i < n; ++i) {
                if (std::isfinite(tThresh[i])) {
                    gr_thA->SetPoint(idx, amps[i], tThresh[i]);
                    gr_thQ->SetPoint(idx, chgs[i], tThresh[i]);
                    idx++;
                }
            }
        }
  
        // 拟合函数（与参考代码一致）
        const string fitForm = "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+"
                         "[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+"
                         "[7]/sqrt(x)+[8]/sqrt(x)/x+[9]/sqrt(x)/x/x";

        // // 拟合函数（与参考代码一致）
        // const string fitForm = "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+"
        //                  "[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+"
        //                  "[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x+[9]*x*x*x*x*x*x*x*x*x";

        // // 拟合函数（与参考代码一致）
        // const string fitForm = "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+"
        //                  "[5]*x*x*x*x*x+[6]*x*x*x*x*x*x";

        // CFD vs Amplitude 拟合
        const double R_l = 1.0;
        const double R_r = 1.0;
        TF1* fitCfdA = new TF1(Form("fF2_cfdA_%dV", voltage), fitForm.c_str(),
                               R_l*max(h2A_lo, 1.0), R_r*h2A_hi);
        // gr_cfdA->Fit(fitCfdA, "RQ");

        // CFD vs Charge 拟合
        TF1* fitCfdQ = new TF1(Form("fF2_cfdQ_%dV", voltage), fitForm.c_str(),
                               R_l*max(h2Q_lo, 100.0), R_r*h2Q_hi);
        // gr_cfdQ->Fit(fitCfdQ, "RQ");

        // Threshold vs Amplitude 拟合
        TF1* fitThA = new TF1(Form("fF2_thA_%dV", voltage), fitForm.c_str(),
                              R_l*max(h2A_lo, 1.0), R_r*h2A_hi);
        // gr_thA->Fit(fitThA, "RQ");

        // Threshold vs Charge 拟合
        TF1* fitThQ = new TF1(Form("fF2_thQ_%dV", voltage), fitForm.c_str(),
                              R_l*max(h2Q_lo, 100.0), R_r*h2Q_hi);
        // gr_thQ->Fit(fitThQ, "RQ");

        // // zyd: 直接对TGraph进行拟合
        // gr_cfdA->Fit(fitCfdA, "RQ");
        // gr_cfdQ->Fit(fitCfdQ, "RQ");
        // gr_thA->Fit(fitThA, "RQ");
        // gr_thQ->Fit(fitThQ, "RQ");

        hCFD_A_ProfileX->Fit(fitCfdA, "RQ");
        hCFD_Q_ProfileX->Fit(fitCfdQ, "RQ");
        hThresh_A_ProfileX->Fit(fitThA, "RQ");
        hThresh_Q_ProfileX->Fit(fitThQ, "RQ");

        // 构建拟合曲线 TGraph 用于绘制
        const int nFitPts = 100;
        TGraph* gFitCfdA = new TGraph(nFitPts);
        TGraph* gFitCfdQ = new TGraph(nFitPts);
        TGraph* gFitThA  = new TGraph(nFitPts);
        TGraph* gFitThQ  = new TGraph(nFitPts);

        for (int i = 0; i < nFitPts; ++i) {
            const double frac = (double)i / (nFitPts - 1);
            const double xA = h2A_lo + frac * (h2A_hi - h2A_lo);
            const double xC = h2Q_lo + frac * (h2Q_hi - h2Q_lo);
            gFitCfdA->SetPoint(i, xA, fitCfdA->Eval(xA));
            gFitCfdQ->SetPoint(i, xC, fitCfdQ->Eval(xC));
            gFitThA->SetPoint(i, xA, fitThA->Eval(xA));
            gFitThQ->SetPoint(i, xC, fitThQ->Eval(xC));
        }

        TCanvas* cFig2 = new TCanvas(Form("cF2_%dV", voltage),
            Form("2D Distributions %dV", voltage), 1600, 1200);
        cFig2->Divide(2, 2);

        // Pad 1: CFD vs Amplitude
        cFig2->cd(1); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        h2_CFD_A->Draw("colz");
        gFitCfdA->SetLineColor(kRed);
        gFitCfdA->SetLineWidth(3);
        gFitCfdA->Draw("L SAME");

        // Pad 2: CFD vs Charge
        cFig2->cd(2); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        h2_CFD_Q->Draw("colz");
        gFitCfdQ->SetLineColor(kRed);
        gFitCfdQ->SetLineWidth(3);
        gFitCfdQ->Draw("L SAME");

        // Pad 3: Threshold vs Amplitude
        cFig2->cd(3); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        h2_Thresh_A->Draw("colz");
        gFitThA->SetLineColor(kRed);
        gFitThA->SetLineWidth(3);
        gFitThA->Draw("L SAME");

        // Pad 4: Threshold vs Charge
        cFig2->cd(4); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        h2_Thresh_Q->Draw("colz");
        gFitThQ->SetLineColor(kRed);
        gFitThQ->SetLineWidth(3);
        gFitThQ->Draw("L SAME");

        cFig2->Write();
        // cFig2->SaveAs(Form("Cali_%dV_Fig2_2D_Distributions.png", voltage));

        // ========== 图3：修正前后时间差分布对比 (2x2) ==========
        // 布局：[1] CFD修正前 | [2] CFD(TA)修正后
        //       [3] 阈值修正前 | [4] 阈值(TA)修正后

        const double caliRange_Tcfd = sigma_num * rmsTcfd;
        const double caliRange_Tthresh = sigma_num * rmsT_thresh;

        TH1F* hT_cfd_before = new TH1F(Form("hF3_tcb_%dV", voltage),
            Form("CFD Before %dV;Time Diff (ns);Counts", voltage), nBinsTcfd, Tcfd_lo, Tcfd_hi);
        TH1F* hT_cfd_TA  = new TH1F(Form("hF3_tca_%dV", voltage),
            Form("CFD After TAC at %dV;Time Diff (ns);Counts", voltage),
            nBinsTcfd, -caliRange_Tcfd, caliRange_Tcfd);
        TH1F* hT_cfd_TQ  = new TH1F(Form("hF3_tcTQ_%dV", voltage),
            Form("CFD After TQC at %dV;Time Diff (ns);Counts", voltage),
            nBinsTcfd, -caliRange_Tcfd, caliRange_Tcfd);


        TH1F* hT_thresh_before = new TH1F(Form("hF3_ttb_%dV", voltage),
            Form("Thresh Before %dV;Time Diff (ns);Counts", voltage), nBinsTthresh, Tthresh_lo, Tthresh_hi);
        TH1F* hT_thresh_TA  = new TH1F(Form("hF3_tta_%dV", voltage),
            Form("Thresh After TAC at %dV;Time Diff (ns);Counts", voltage),
            nBinsTthresh, -caliRange_Tthresh, caliRange_Tthresh);
        TH1F* hT_thresh_TQ  = new TH1F(Form("hF3_ttTQ_%dV", voltage),
            Form("Thresh After TQC at %dV;Time Diff (ns);Counts", voltage),
            nBinsTthresh, -caliRange_Tthresh, caliRange_Tthresh);

        for (int i = 0; i < n; ++i) {
            hT_cfd_before->Fill(tCFD[i]);
            hT_cfd_TA->Fill(tCFD[i] - fitCfdA->Eval(amps[i]));
            hT_cfd_TQ->Fill(tCFD[i] - fitCfdQ->Eval(chgs[i]));
            if (std::isfinite(tThresh[i])) {
                hT_thresh_before->Fill(tThresh[i]);
                hT_thresh_TA->Fill(tThresh[i] - fitThA->Eval(amps[i]));
                hT_thresh_TQ->Fill(tThresh[i] - fitThQ->Eval(chgs[i]));
            }
        }

        TCanvas* cFig3 = new TCanvas(Form("cF3_%dV", voltage),
            Form("Correction Comparison %dV", voltage), 1600, 1200);
        cFig3->Divide(3, 2);

        // Pad 1: CFD 修正前
        cFig3->cd(1); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        hT_cfd_before->Draw();
        {
            double m, s, me, se;
            iterativeGaussianFit(hT_cfd_before, 3, m, s, me, se);
            cfd_sigma_before[voltage] = s;
            cfd_sigmaErr_before[voltage] = se;
        }

        // Pad 2: CFD TA修正后
        cFig3->cd(2); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        hT_cfd_TA->Draw();
        {
            double m, s, me, se;
            iterativeGaussianFit(hT_cfd_TA, 3, m, s, me, se);
            cfd_sigma_afterTA[voltage] = s;
            cfd_sigmaErr_afterTA[voltage] = se;
        }

        // Pad 3: CFD TQ修正后
        cFig3->cd(3); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        hT_cfd_TQ->Draw();
        {
            double m, s, me, se;
            iterativeGaussianFit(hT_cfd_TQ, 3, m, s, me, se);
            cfd_sigma_afterTQ[voltage] = s;
            cfd_sigmaErr_afterTQ[voltage] = se;
        }

        // Pad 4: 阈值 修正前
        cFig3->cd(4); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        hT_thresh_before->Draw();
        {
            double m, s, me, se;
            iterativeGaussianFit(hT_thresh_before, 3, m, s, me, se);
            thresh_sigma_before[voltage] = s;
            thresh_sigmaErr_before[voltage] = se;
        }

        // Pad 5: 阈值 TA修正后
        cFig3->cd(5); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        hT_thresh_TA->Draw();
        {
            double m, s, me, se;
            iterativeGaussianFit(hT_thresh_TA, 3, m, s, me, se);
            thresh_sigma_afterTA[voltage] = s;
            thresh_sigmaErr_afterTA[voltage] = se;
        }
        // Pad 6: 阈值 TQ修正后
        cFig3->cd(6); gPad->SetMargin(0.12, 0.05, 0.12, 0.05);
        hT_thresh_TQ->Draw();
        {
            double m, s, me, se;
            iterativeGaussianFit(hT_thresh_TQ, 3, m, s, me, se);
            thresh_sigma_afterTQ[voltage] = s;
            thresh_sigmaErr_afterTQ[voltage] = se;
        }
        

        cFig3->Write();
        // cFig3->SaveAs(Form("Cali_%dV_Fig3_Correction_Comparison.png", voltage));

        // ===== 计算 TQ 修正后的时间分辨（不画图，仅用于汇总） =====
        {
            TH1F* hTC_tq = new TH1F(Form("hF4_tcTQ_%dV", voltage), "", nBinsTcfd, -caliRange_Tcfd, caliRange_Tcfd);
            TH1F* hTT_tq = new TH1F(Form("hF4_ttTQ_%dV", voltage), "", nBinsTthresh, -caliRange_Tthresh, caliRange_Tthresh);
            for (int i = 0; i < n; ++i) {
                hTC_tq->Fill(tCFD[i] - fitCfdQ->Eval(chgs[i]));
                if (std::isfinite(tThresh[i])) {
                    hTT_tq->Fill(tThresh[i] - fitThQ->Eval(chgs[i]));
                }
            }
            
            double m, s, me, se;
            iterativeGaussianFit(hTC_tq, 3, m, s, me, se);
            cfd_sigma_afterTQ[voltage] = s;
            cfd_sigmaErr_afterTQ[voltage] = se;
            iterativeGaussianFit(hTT_tq, 3, m, s, me, se);
            thresh_sigma_afterTQ[voltage] = s;
            thresh_sigmaErr_afterTQ[voltage] = se;
            delete hTC_tq;
            delete hTT_tq;
        }

        cout << "  CFD:  before=" << cfd_sigma_before[voltage]
             << "  TA=" << cfd_sigma_afterTA[voltage]
             << "  TQ=" << cfd_sigma_afterTQ[voltage] << " ns" << endl;
        cout << "  Thresh: before=" << thresh_sigma_before[voltage]
             << "  TA=" << thresh_sigma_afterTA[voltage]
             << "  TQ=" << thresh_sigma_afterTQ[voltage] << " ns" << endl;

        // 释放图2相关对象内存 - 修复内存泄漏
        delete gr_cfdA;
        delete gr_cfdQ;
        delete gr_thA;
        delete gr_thQ;
        delete fitCfdA;
        delete fitCfdQ;
        delete fitThA;
        delete fitThQ;
        delete gFitCfdA;
        delete gFitCfdQ;
        delete gFitThA;
        delete gFitThQ;
        delete h2_CFD_A;
        delete h2_CFD_Q;
        delete h2_Thresh_A;
        delete h2_Thresh_Q;

        delete hCFD_A_ProfileX;
        delete hCFD_Q_ProfileX;
        delete hThresh_A_ProfileX;
        delete hThresh_Q_ProfileX;

        // 返回根目录
        outputFile->cd();
    }

   // ========== 分辨率 vs 电压图（一级目录） ==========
    outputFile->cd();
    const int nv = static_cast<int>(voltages.size());

    // CFD 定时
    TGraphErrors* gr_cfd_before = new TGraphErrors(nv);
    TGraphErrors* gr_cfd_afterTA = new TGraphErrors(nv);
    TGraphErrors* gr_cfd_afterTQ = new TGraphErrors(nv);

    // 阈值定时
    TGraphErrors* gr_th_before = new TGraphErrors(nv);
    TGraphErrors* gr_th_afterTA = new TGraphErrors(nv);
    TGraphErrors* gr_th_afterTQ = new TGraphErrors(nv);

    // 用于计算全局坐标轴范围
    double x_min =  1e30, x_max = -1e30;
    double y_min =  1e30, y_max = -1e30;

    for (int i = 0; i < nv; ++i) {
        const int v = voltages[i];
        const double x = (double)v;

        // 纳秒 → 皮秒：乘以 1000
        double cfd_b    = (cfd_sigma_before.count(v)    ? cfd_sigma_before[v]    : 0) * 1000.0;
        double cfd_bErr = (cfd_sigmaErr_before.count(v)  ? cfd_sigmaErr_before[v]  : 0) * 1000.0;
        double cfd_aTA    = (cfd_sigma_afterTA.count(v)    ? cfd_sigma_afterTA[v]    : 0) * 1000.0;
        double cfd_aTAErr = (cfd_sigmaErr_afterTA.count(v)  ? cfd_sigmaErr_afterTA[v]  : 0) * 1000.0;
        double cfd_aTQ    = (cfd_sigma_afterTQ.count(v)    ? cfd_sigma_afterTQ[v]    : 0) * 1000.0;
        double cfd_aTQErr = (cfd_sigmaErr_afterTQ.count(v)  ? cfd_sigmaErr_afterTQ[v]  : 0) * 1000.0;

        double th_b    = (thresh_sigma_before.count(v)    ? thresh_sigma_before[v]    : 0) * 1000.0;
        double th_bErr = (thresh_sigmaErr_before.count(v)  ? thresh_sigmaErr_before[v]  : 0) * 1000.0;
        double th_aTA    = (thresh_sigma_afterTA.count(v)    ? thresh_sigma_afterTA[v]    : 0) * 1000.0;
        double th_aTAErr = (thresh_sigmaErr_afterTA.count(v)  ? thresh_sigmaErr_afterTA[v]  : 0) * 1000.0;
        double th_aTQ    = (thresh_sigma_afterTQ.count(v)    ? thresh_sigma_afterTQ[v]    : 0) * 1000.0;
        double th_aTQErr = (thresh_sigmaErr_afterTQ.count(v)  ? thresh_sigmaErr_afterTQ[v]  : 0) * 1000.0;

        gr_cfd_before->SetPoint(i, x, cfd_b);
        gr_cfd_before->SetPointError(i, 0, cfd_bErr);
        gr_cfd_afterTA->SetPoint(i, x, cfd_aTA);
        gr_cfd_afterTA->SetPointError(i, 0, cfd_aTAErr);
        gr_cfd_afterTQ->SetPoint(i, x, cfd_aTQ);
        gr_cfd_afterTQ->SetPointError(i, 0, cfd_aTQErr);

        gr_th_before->SetPoint(i, x, th_b);
        gr_th_before->SetPointError(i, 0, th_bErr);
        gr_th_afterTA->SetPoint(i, x, th_aTA);
        gr_th_afterTA->SetPointError(i, 0, th_aTAErr);
        gr_th_afterTQ->SetPoint(i, x, th_aTQ);
        gr_th_afterTQ->SetPointError(i, 0, th_aTQErr);

        // 统计 x 范围
        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;

        // 统计 y 范围（含误差棒）
        auto updateY = [&](double val, double err) {
            if (val - err < y_min) y_min = val - err;
            if (val + err > y_max) y_max = val + err;
        };
        updateY(cfd_b,    cfd_bErr);
        updateY(cfd_aTA,  cfd_aTAErr);
        updateY(cfd_aTQ,  cfd_aTQErr);
        updateY(th_b,     th_bErr);
        updateY(th_aTA,   th_aTAErr);
        updateY(th_aTQ,   th_aTQErr);
    }

    // 计算坐标轴范围：0.9倍最小值 ~ 1.1倍最大值
    double X_MIN = x_min * 0.95;
    double X_MAX = x_max * 1.1;
    double Y_MIN = y_min * 0.8;
    double Y_MAX = y_max * 1.4;
    if (y_min <= 0) Y_MIN = 0;

    // ========== 设置样式 ==========
    // CFD：空心标记
    gr_cfd_before->SetMarkerStyle(24);   gr_cfd_before->SetMarkerColor(kBlue);
    gr_cfd_before->SetLineColor(kBlue);  gr_cfd_before->SetLineWidth(3);
    gr_cfd_afterTA->SetMarkerStyle(25);  gr_cfd_afterTA->SetMarkerColor(kRed);
    gr_cfd_afterTA->SetLineColor(kRed);  gr_cfd_afterTA->SetLineWidth(3);
    gr_cfd_afterTQ->SetMarkerStyle(26);  gr_cfd_afterTQ->SetMarkerColor(kGreen+2);
    gr_cfd_afterTQ->SetLineColor(kGreen+2); gr_cfd_afterTQ->SetLineWidth(3);

    // 阈值：实心标记 + 虚线
    gr_th_before->SetMarkerStyle(20);    gr_th_before->SetMarkerColor(kBlue);
    gr_th_before->SetLineColor(kBlue);   gr_th_before->SetLineWidth(3);
    gr_th_before->SetLineStyle(2);
    gr_th_afterTA->SetMarkerStyle(21);   gr_th_afterTA->SetMarkerColor(kRed);
    gr_th_afterTA->SetLineColor(kRed);   gr_th_afterTA->SetLineWidth(3);
    gr_th_afterTA->SetLineStyle(2);
    gr_th_afterTQ->SetMarkerStyle(22);   gr_th_afterTQ->SetMarkerColor(kGreen+2);
    gr_th_afterTQ->SetLineColor(kGreen+2); gr_th_afterTQ->SetLineWidth(3);
    gr_th_afterTQ->SetLineStyle(2);

    // ========== 图1：CFD 时间分辨 vs 电压 ==========
    TCanvas* c_res_cfd = new TCanvas("c_res_cfd", "CFD Time Resolution vs Voltage", 1000, 700);
    c_res_cfd->SetGrid();
    TMultiGraph* mg_cfd = new TMultiGraph();
    mg_cfd->SetTitle("CFD Time Resolution vs HV;Voltage (V);Time Resolution (ps)");
    mg_cfd->Add(gr_cfd_before, "LP");
    mg_cfd->Add(gr_cfd_afterTA, "LP");
    mg_cfd->Add(gr_cfd_afterTQ, "LP");
    mg_cfd->Draw("A");
    mg_cfd->GetXaxis()->SetLimits(X_MIN, X_MAX);
    mg_cfd->GetXaxis()->SetRangeUser(X_MIN, X_MAX);
    mg_cfd->GetYaxis()->SetLimits(Y_MIN, Y_MAX);
    mg_cfd->GetYaxis()->SetRangeUser(Y_MIN, Y_MAX);
    mg_cfd->GetXaxis()->SetTitleOffset(1.1);
    mg_cfd->GetYaxis()->SetTitleOffset(1.1);
    c_res_cfd->Modified();
    c_res_cfd->Update();

    TLegend* leg_cfd = new TLegend(0.55, 0.7, 0.9, 0.9);
    leg_cfd->SetTextSize(0.03);
    leg_cfd->AddEntry(gr_cfd_before, "CFD Before Correction", "lp");
    leg_cfd->AddEntry(gr_cfd_afterTA, "CFD After TA Correction", "lp");
    leg_cfd->AddEntry(gr_cfd_afterTQ, "CFD After TQ Correction", "lp");
    leg_cfd->Draw();
    c_res_cfd->Write();
    // c_res_cfd->SaveAs(Form("%s_CFD_TR_vs_HV.png", outputFileName));

    // ========== 图2：阈值定时 时间分辨 vs 电压 ==========
    TCanvas* c_res_th = new TCanvas("c_res_th", "Threshold Time Resolution vs Voltage", 1000, 700);
    c_res_th->SetGrid();
    TMultiGraph* mg_th = new TMultiGraph();
    mg_th->SetTitle("Threshold Time Resolution vs HV;Voltage (V);Time Resolution (ps)");
    mg_th->Add(gr_th_before, "LP");
    mg_th->Add(gr_th_afterTA, "LP");
    mg_th->Add(gr_th_afterTQ, "LP");
    mg_th->Draw("A");
    mg_th->GetXaxis()->SetLimits(X_MIN, X_MAX);
    mg_th->GetXaxis()->SetRangeUser(X_MIN, X_MAX);
    mg_th->GetYaxis()->SetLimits(Y_MIN, Y_MAX);
    mg_th->GetYaxis()->SetRangeUser(Y_MIN, Y_MAX);
    mg_th->GetXaxis()->SetTitleOffset(1.1);
    mg_th->GetYaxis()->SetTitleOffset(1.1);
    c_res_th->Modified();
    c_res_th->Update();

    TLegend* leg_th = new TLegend(0.55, 0.7, 0.9, 0.9);
    leg_th->SetTextSize(0.03);
    leg_th->AddEntry(gr_th_before, "Threshold Before Correction", "lp");
    leg_th->AddEntry(gr_th_afterTA, "Threshold After TA Correction", "lp");
    leg_th->AddEntry(gr_th_afterTQ, "Threshold After TQ Correction", "lp");
    leg_th->Draw();
    c_res_th->Write();
    // c_res_th->SaveAs(Form("%s_Threshold_TR_vs_HV.png", outputFileName));

    // ========== 图3：CFD + 阈值 六条曲线合并 ==========
    TCanvas* c_res_combined = new TCanvas("c_res_combined",
        "Combined Time Resolution vs Voltage", 1000, 700);
    c_res_combined->SetGrid();
    TMultiGraph* mg_combined = new TMultiGraph();
    mg_combined->SetTitle("Time Resolution vs HV (CFD & Threshold);Voltage (V);Time Resolution (ps)");
    mg_combined->Add(gr_cfd_before, "LP");
    mg_combined->Add(gr_cfd_afterTA, "LP");
    mg_combined->Add(gr_cfd_afterTQ, "LP");
    mg_combined->Add(gr_th_before, "LP");
    mg_combined->Add(gr_th_afterTA, "LP");
    mg_combined->Add(gr_th_afterTQ, "LP");
    mg_combined->Draw("A");
    mg_combined->GetXaxis()->SetLimits(X_MIN, X_MAX);
    mg_combined->GetXaxis()->SetRangeUser(X_MIN, X_MAX);
    mg_combined->GetYaxis()->SetLimits(Y_MIN, Y_MAX);
    mg_combined->GetYaxis()->SetRangeUser(Y_MIN, Y_MAX);
    mg_combined->GetXaxis()->SetTitleOffset(1.1);
    mg_combined->GetYaxis()->SetTitleOffset(1.1);
    c_res_combined->Modified();
    c_res_combined->Update();

    TLegend* leg_combined = new TLegend(0.45, 0.65, 0.9, 0.9);
    leg_combined->SetTextSize(0.025);
    leg_combined->AddEntry(gr_cfd_before,  "CFD - Before Correction",       "lp");
    leg_combined->AddEntry(gr_cfd_afterTA, "CFD - After TA Correction",     "lp");
    leg_combined->AddEntry(gr_cfd_afterTQ, "CFD - After TQ Correction",     "lp");
    leg_combined->AddEntry(gr_th_before,   "Threshold - Before Correction",  "lp");
    leg_combined->AddEntry(gr_th_afterTA,  "Threshold - After TA Correction","lp");
    leg_combined->AddEntry(gr_th_afterTQ,  "Threshold - After TQ Correction","lp");
    leg_combined->Draw();
    c_res_combined->Write();
    // c_res_combined->SaveAs(Form("%s_Combined_TR_vs_HV.png", outputFileName));


    //修改这段代码，增加第三幅度，把六条曲线放在同一幅图中进行对比，用合适的图例和颜色区分不同曲线，并保存； 作图时要求横、纵坐标范围设置为：最小值为数据最小值的0.9倍，最大值为数据最大值的1.1倍 ，设置范围参考下面的代码


    // // 输出汇总表
    // cout << "\n\n========== 时间分辨汇总 ==========" << endl;
    // cout << "电压(V)\tCFD_before(ns)\tCFD_TA(ns)\tCFD_TQ(ns)\tThresh_before(ns)\tThresh_TA(ns)\tThresh_TQ(ns)" << endl;
    // for (int v : voltages) {
    //     cout << v << "\t"
    //          << (cfd_sigma_before.count(v) ? cfd_sigma_before[v] : 0) << " +/- "
    //          << (cfd_sigmaErr_before.count(v) ? cfd_sigmaErr_before[v] : 0) << "\t"
    //          << (cfd_sigma_afterTA.count(v) ? cfd_sigma_afterTA[v] : 0) << " +/- "
    //          << (cfd_sigmaErr_afterTA.count(v) ? cfd_sigmaErr_afterTA[v] : 0) << "\t"
    //          << (cfd_sigma_afterTQ.count(v) ? cfd_sigma_afterTQ[v] : 0) << " +/- "
    //          << (cfd_sigmaErr_afterTQ.count(v) ? cfd_sigmaErr_afterTQ[v] : 0) << "\t"
    //          << (thresh_sigma_before.count(v) ? thresh_sigma_before[v] : 0) << " +/- "
    //          << (thresh_sigmaErr_before.count(v) ? thresh_sigmaErr_before[v] : 0) << "\t"
    //          << (thresh_sigma_afterTA.count(v) ? thresh_sigma_afterTA[v] : 0) << " +/- "
    //          << (thresh_sigmaErr_afterTA.count(v) ? thresh_sigmaErr_afterTA[v] : 0) << "\t"
    //          << (thresh_sigma_afterTQ.count(v) ? thresh_sigma_afterTQ[v] : 0) << " +/- "
    //          << (thresh_sigmaErr_afterTQ.count(v) ? thresh_sigmaErr_afterTQ[v] : 0)
    //          << endl;
    // }
}



// 优化后的二维图创建 - 使用查找表替代if-else链
void create2DPlotsOptimized(const vector<int>& voltages,
                           const map<int, FitParameters>& FitParamsMap,
                           const map<int, vector<double>>& amplitudes,
                           const map<int, vector<double>>& rise_times,
                           const map<int, vector<double>>& fwhms,
                           const map<int, vector<double>>& charges,
                           const map<int, vector<double>>& timing_diffs_const,
                           TFile* outputFile) {

    outputFile->cd();

    vector<pair<string, pair<string, string>>> plot_combinations = {
        {"Amplitude_vs_RiseTime", {"Amplitude (mV)", "Rise Time (ns)"}},
        {"Amplitude_vs_FWHM", {"Amplitude (mV)", "FWHM (ns)"}},
        {"Amplitude_vs_Charge", {"Amplitude (mV)", "Charge (Q_{e})"}},
        {"Amplitude_vs_Timing", {"Amplitude (mV)", "Time Difference (ns)"}},
        {"RiseTime_vs_FWHM", {"Rise Time (ns)", "FWHM (ns)"}},
        {"RiseTime_vs_Charge", {"Rise Time (ns)", "Charge (Q_{e})"}},
        {"RiseTime_vs_Timing", {"Rise Time (ns)", "Time Difference (ns)"}},
        {"FWHM_vs_Charge", {"FWHM (ns)", "Charge (Q_{e})"}},
        {"FWHM_vs_Timing", {"FWHM (ns)", "Time Difference (ns)"}},
        {"Charge_vs_Timing", {"Charge (Q_{e})", "Time Difference (ns)"}}
    };

    for (const auto& combo : plot_combinations) {
        const string& plot_name = combo.first;
        const string& x_title = combo.second.first;
        const string& y_title = combo.second.second;

        TDirectory* combo_dir = outputFile->mkdir(plot_name.c_str());
        combo_dir->cd();

        map<int, TH2F*> h2D;

        for (int voltage : voltages) {
            const auto& params = FitParamsMap.at(voltage);

            const double Num_Sigma = 10.0;
            double x_min, x_max, y_min, y_max;

            // 优化：使用结构化绑定和映射表减少重复代码
            // 根据plot_name获取对应的参数
            double mean_x, sigma_x, mean_y, sigma_y;

            if (plot_name == "Amplitude_vs_RiseTime") {
                mean_x = params.mean_amp; sigma_x = params.sigma_amp;
                mean_y = params.mean_rt; sigma_y = params.sigma_rt;
            } else if (plot_name == "Amplitude_vs_FWHM") {
                mean_x = params.mean_amp; sigma_x = params.sigma_amp;
                mean_y = params.mean_fwhm; sigma_y = params.sigma_fwhm;
            } else if (plot_name == "Amplitude_vs_Charge") {
                mean_x = params.mean_amp; sigma_x = params.sigma_amp;
                mean_y = params.mean_charge; sigma_y = params.sigma_charge;
            } else if (plot_name == "Amplitude_vs_Timing") {
                mean_x = params.mean_amp; sigma_x = params.sigma_amp;
                mean_y = params.mean_timing_CFD; sigma_y = params.sigma_timing_CFD;
            } else if (plot_name == "RiseTime_vs_FWHM") {
                mean_x = params.mean_rt; sigma_x = params.sigma_rt;
                mean_y = params.mean_fwhm; sigma_y = params.sigma_fwhm;
            } else if (plot_name == "RiseTime_vs_Charge") {
                mean_x = params.mean_rt; sigma_x = params.sigma_rt;
                mean_y = params.mean_charge; sigma_y = params.sigma_charge;
            } else if (plot_name == "RiseTime_vs_Timing") {
                mean_x = params.mean_rt; sigma_x = params.sigma_rt;
                mean_y = params.mean_timing_CFD; sigma_y = params.sigma_timing_CFD;
            } else if (plot_name == "FWHM_vs_Charge") {
                mean_x = params.mean_fwhm; sigma_x = params.sigma_fwhm;
                mean_y = params.mean_charge; sigma_y = params.sigma_charge;
            } else if (plot_name == "FWHM_vs_Timing") {
                mean_x = params.mean_fwhm; sigma_x = params.sigma_fwhm;
                mean_y = params.mean_timing_CFD; sigma_y = params.sigma_timing_CFD;
            } else { // Charge_vs_Timing
                mean_x = params.mean_charge; sigma_x = params.sigma_charge;
                mean_y = params.mean_timing_CFD; sigma_y = params.sigma_timing_CFD;
            }

            x_min = max(0.0, mean_x - Num_Sigma * sigma_x);
            x_max = mean_x + Num_Sigma * sigma_x;
            y_min = max(0.0, mean_y - Num_Sigma * sigma_y);
            y_max = mean_y + Num_Sigma * sigma_y;

            // x_min = mean_x - Num_Sigma * sigma_x;
            // x_max = mean_x + Num_Sigma * sigma_x;
            // y_min = mean_y - Num_Sigma * sigma_y;
            // y_max = mean_y + Num_Sigma * sigma_y;

            if (plot_name.find("Timing") != std::string::npos)
            {
                y_min = mean_y - Num_Sigma * sigma_y;
            }

            // if ((plot_name.find("Charge") != std::string::npos) || (plot_name.find("Amplitude") != std::string::npos))
            // {
            //     x_max = mean_x + 5 * sigma_x;
            //     y_max = mean_y + 5 * sigma_y;
            // }
            
            

            const double x_range = x_max - x_min;
            const double y_range = y_max - y_min;

            const double binX_width = (x_range < 20) ? 0.01 : (x_range < 1e4) ? 2 : 1e4;
            const double binY_width = (y_range < 20) ? 0.01 : (y_range < 1e4) ? 2 : 1e4;

            const int nBinsX = max(static_cast<int>(ceil((x_max - x_min) / binX_width)), 50);
            const int nBinsY = max(static_cast<int>(ceil((y_max - y_min) / binY_width)), 50);

            h2D[voltage] = new TH2F(Form("h2D_%s_%dV", plot_name.c_str(), voltage),
                                   Form("%s at %dV;%s;%s", plot_name.c_str(), voltage,
                                        x_title.c_str(), y_title.c_str()),
                                   nBinsX, x_min, x_max, nBinsY, y_min, y_max);

            const auto& amps = amplitudes.at(voltage);
            const auto& rts = rise_times.at(voltage);
            const auto& fwms = fwhms.at(voltage);
            const auto& chgs = charges.at(voltage);
            const auto& tdiffs = timing_diffs_const.at(voltage);

            const size_t data_size = amps.size();
            for (size_t j = 0; j < data_size; ++j) {
                double x, y;

                if (plot_name == "Amplitude_vs_RiseTime") {
                    x = amps[j]; y = rts[j];
                } else if (plot_name == "Amplitude_vs_FWHM") {
                    x = amps[j]; y = fwms[j];
                } else if (plot_name == "Amplitude_vs_Charge") {
                    x = amps[j]; y = chgs[j];
                } else if (plot_name == "Amplitude_vs_Timing") {
                    x = amps[j]; y = tdiffs[j];
                } else if (plot_name == "RiseTime_vs_FWHM") {
                    x = rts[j]; y = fwms[j];
                } else if (plot_name == "RiseTime_vs_Charge") {
                    x = rts[j]; y = chgs[j];
                } else if (plot_name == "RiseTime_vs_Timing") {
                    x = rts[j]; y = tdiffs[j];
                } else if (plot_name == "FWHM_vs_Charge") {
                    x = fwms[j]; y = chgs[j];
                } else if (plot_name == "FWHM_vs_Timing") {
                    x = fwms[j]; y = tdiffs[j];
                } else if (plot_name == "Charge_vs_Timing") {
                    x = chgs[j]; y = tdiffs[j];
                }

                h2D[voltage]->Fill(x, y);
            }

            h2D[voltage]->Write();
        }

        TCanvas* canvas = new TCanvas(Form("c_%s", plot_name.c_str()),
                                     Form("%s", plot_name.c_str()), 1600, 1600);
        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

        bool first = true;
        for (int voltage : voltages) {
            if (first) {
                h2D[voltage]->Draw("COLZ");
                first = false;
            } else {
                h2D[voltage]->Draw("COLZ SAME");
            }
            legend->AddEntry(h2D[voltage], Form("%dV", voltage), "f");
        }

        legend->Draw();
        canvas->Write();

        delete canvas;
        delete legend;

        for (int voltage : voltages) {
            delete h2D[voltage];
        }

        outputFile->cd();
    }
}

// ========== 修改后的主分析函数 ==========
void analyzeDataSerial() {
    // gStyle->SetOptStat(111111);
    // gStyle->SetOptFit(1111);
    // gStyle->SetPalette(55);
    // gStyle->SetOptStat(111111);
    gStyle->SetOptFit(111);
    // gStyle->SetPalette(55);

    // zyd
    // vector<int> voltages = {1000};
    // vector<int> voltages = { 1900, 1925, 1950, 1975, 2000, 2025, 2050, 2075, 2100, 2125, 2150, 2175, 2200, 2225, 2250, 2275, 2300, 2325, 2350, 2375, 2400, 2450, 2500, 2550};

    vector<int> voltages = GetVoltages(csvFolder);
    // // 打印结果
    // for (auto v : voltages) {
    //     cout << v << "V" << endl;
    // }

    // TFile* outputFile = new TFile("RPC_PMT_Analysis.root", "RECREATE"); //outputFileName
    TFile* outputFile = new TFile(Form("%s.root", outputFileName), "RECREATE");

    map<int, vector<double>> amplitudes, rise_times, fwhms, charges;
    map<int, vector<double>> timing_diffs_const, timing_diffs_thresh;
    map<int, FitParameters> fit_params_map;

    TDirectory* dir_amp = outputFile->mkdir("amplitude");
    TDirectory* dir_rt = outputFile->mkdir("rise_time");
    TDirectory* dir_fwhm = outputFile->mkdir("fwhm");
    TDirectory* dir_charge = outputFile->mkdir("charge");
    TDirectory* dir_timing_const = outputFile->mkdir("timing_const");
    TDirectory* dir_timing_thresh = outputFile->mkdir("timing_thresh");
    // // 创建FitParameters目录
    // TDirectory* dir_fit = outputFile->mkdir("FitParameters");

    // 创建TTree存储拟合参数
    TTree* fit_tree = new TTree("tree_fit_parameters", "Fit Parameters for Each Voltage");
    int voltage_var;
    double mean_amp, sigma_amp, mean_rt, sigma_rt, mean_fwhm, sigma_fwhm;
    double mean_charge, sigma_charge, mean_timing_CFD, sigma_timing_CFD;
    double mean_timing_thresh, sigma_timing_thresh;

    fit_tree->Branch("voltage", &voltage_var, "voltage/I");
    fit_tree->Branch("mean_amp", &mean_amp, "mean_amp/D");
    fit_tree->Branch("sigma_amp", &sigma_amp, "sigma_amp/D");
    fit_tree->Branch("mean_rt", &mean_rt, "mean_rt/D");
    fit_tree->Branch("sigma_rt", &sigma_rt, "sigma_rt/D");
    fit_tree->Branch("mean_fwhm", &mean_fwhm, "mean_fwhm/D");
    fit_tree->Branch("sigma_fwhm", &sigma_fwhm, "sigma_fwhm/D");
    fit_tree->Branch("mean_charge", &mean_charge, "mean_charge/D");
    fit_tree->Branch("sigma_charge", &sigma_charge, "sigma_charge/D");
    fit_tree->Branch("mean_timing_CFD", &mean_timing_CFD, "mean_timing_CFD/D");
    fit_tree->Branch("sigma_timing_CFD", &sigma_timing_CFD, "sigma_timing_CFD/D");
    fit_tree->Branch("mean_timing_thresh", &mean_timing_thresh, "mean_timing_thresh/D");
    fit_tree->Branch("sigma_timing_thresh", &sigma_timing_thresh, "sigma_timing_thresh/D");

    // 创建TTree存储TAx修正后的时间分辨率
    TTree* TAC_tree = new TTree("tree_fitTAC_parameters", "Fit Parameters After TAC for Each Voltage");

    double sigma_cfd_before, sigmaErr_cfd_before;
    double sigma_cfd_afterTA, sigmaErr_cfd_afterTA;
    double sigma_cfd_afterTQ, sigmaErr_cfd_afterTQ;

    double sigma_thresh_before, sigmaErr_thresh_before;
    double sigma_thresh_afterTA, sigmaErr_thresh_afterTA;
    double sigma_thresh_afterTQ, sigmaErr_thresh_afterTQ;

    TAC_tree->Branch("sigma_cfd_before", &sigma_cfd_before, "sigma_cfd_before/D");
    TAC_tree->Branch("sigmaErr_cfd_before", &sigmaErr_cfd_before, "sigmaErr_cfd_before/D");
    TAC_tree->Branch("sigma_cfd_afterTA", &sigma_cfd_afterTA, "sigma_cfd_afterTA/D");
    TAC_tree->Branch("sigmaErr_cfd_afterTA", &sigmaErr_cfd_afterTA, "sigmaErr_cfd_afterTA/D");
    TAC_tree->Branch("sigma_cfd_afterTQ", &sigma_cfd_afterTQ, "sigma_cfd_afterTQ/D");
    TAC_tree->Branch("sigmaErr_cfd_afterTQ", &sigmaErr_cfd_afterTQ, "sigmaErr_cfd_afterTQ/D");

    TAC_tree->Branch("sigma_thresh_before", &sigma_thresh_before, "sigma_thresh_before/D");
    TAC_tree->Branch("sigmaErr_thresh_before", &sigmaErr_thresh_before, "sigmaErr_thresh_before/D");
    TAC_tree->Branch("sigma_thresh_afterTA", &sigma_thresh_afterTA, "sigma_thresh_afterTA/D");
    TAC_tree->Branch("sigmaErr_thresh_afterTA", &sigmaErr_thresh_afterTA, "sigmaErr_thresh_afterTA/D");
    TAC_tree->Branch("sigma_thresh_afterTQ", &sigma_thresh_afterTQ, "sigma_thresh_afterTQ/D");
    TAC_tree->Branch("sigmaErr_thresh_afterTQ", &sigmaErr_thresh_afterTQ, "sigmaErr_thresh_afterTQ/D");

    int total_files = 0;
    for (int voltage : voltages) {
        string folder_name = to_string(voltage);

        TTree* tree = new TTree(Form("tree_%dV", voltage),
                               Form("RPC/PMT data at %dV", voltage));

        cout << "处理电压: " << voltage << " V" << endl;

        processVoltageSerial(voltage, folder_name,
                           amplitudes[voltage],
                           rise_times[voltage],
                           fwhms[voltage],
                           charges[voltage],
                           timing_diffs_const[voltage],
                           timing_diffs_thresh[voltage],
                           tree);

        total_files += amplitudes[voltage].size();

        // 第一次创建直方图并进行高斯拟合
        // zyd
        const double time_bin_width = 0.008;

        const double RT_min = 0, RT_max = 3;
        const double fwhm_min = 0.5, fwhm_max = 4;
        const double DT_min = -5, DT_max = 10;

        const double A_bin_width = 1.5; // 幅值的bin宽度设置为2 mV
        const double Q_bin_width = 2e4; // 电荷的bin宽度设置为1e3 e
        const double hA_min = 0, hA_max = 1500; // 幅值范围设置为0-1500 mV
        const double hQ_min = 0, hQ_max = 3e7; // 电荷范围设置为0-3e7 e

        const int RT_bins = static_cast<int>((RT_max - RT_min) / (0.003)); // Rise Time的bin数量设置为原来的2倍
        const int fwhm_bins = static_cast<int>((fwhm_max - fwhm_min) / time_bin_width);
        const int DT_bins = static_cast<int>((DT_max - DT_min) / time_bin_width);

        const int A_bins = static_cast<int>((hA_max - hA_min) / A_bin_width);
        const int Q_bins = static_cast<int>((hQ_max - hQ_min) / Q_bin_width);

        TH1F* h_amplitude = new TH1F(Form("h_amp_%dV", voltage),
                                    Form("Amplitude at %dV;Amplitude (mV);Counts", voltage),
                                    A_bins, hA_min, hA_max);
        TH1F* h_rise_time = new TH1F(Form("h_rt_%dV", voltage),
                                    Form("Rise Time at %dV;Rise Time (ns);Counts", voltage),
                                    RT_bins, RT_min, RT_max);
        TH1F* h_fwhm = new TH1F(Form("h_fwhm_%dV", voltage),
                               Form("FWHM at %dV;FWHM (ns);Counts", voltage),
                               fwhm_bins, fwhm_min, fwhm_max);
        TH1F* h_charge = new TH1F(Form("h_charge_%dV", voltage),
                                 Form("Charge at %dV;Charge (e);Counts", voltage),
                                 Q_bins, hQ_min, hQ_max);
        TH1F* h_timing_const = new TH1F(Form("h_timing_CFD_%dV", voltage),
                                       Form("Timing (CF=%.2f) at %dV;Time Difference (ns);Counts", RPC_CFD_FRACTION, voltage),
                                       DT_bins, DT_min, DT_max);
        TH1F* h_timing_thresh = new TH1F(Form("h_timing_thresh_%dV", voltage),
                                        Form("Timing (Threshold=%.1f mV) at %dV;Time Difference (ns);Counts", RPC_THRESHOLD, voltage),
                                        DT_bins, DT_min, DT_max);

        const size_t amp_size = amplitudes[voltage].size();
        for (size_t i = 0; i < amp_size; ++i) {
            h_amplitude->Fill(amplitudes[voltage][i]);
            h_rise_time->Fill(rise_times[voltage][i]);
            h_fwhm->Fill(fwhms[voltage][i]);
            h_charge->Fill(charges[voltage][i]);
            h_timing_const->Fill(timing_diffs_const[voltage][i]);

            if (std::isfinite(timing_diffs_thresh[voltage][i])) {
                h_timing_thresh->Fill(timing_diffs_thresh[voltage][i]);
            }
        }

        fit_params_map[voltage] = FitParameters();
        gaussianFit(h_amplitude, fit_params_map[voltage].mean_amp, fit_params_map[voltage].sigma_amp);
        gaussianFit(h_rise_time, fit_params_map[voltage].mean_rt, fit_params_map[voltage].sigma_rt);
        gaussianFit(h_fwhm, fit_params_map[voltage].mean_fwhm, fit_params_map[voltage].sigma_fwhm);
        gaussianFit(h_charge, fit_params_map[voltage].mean_charge, fit_params_map[voltage].sigma_charge);
        gaussianFit(h_timing_const, fit_params_map[voltage].mean_timing_CFD, fit_params_map[voltage].sigma_timing_CFD);
        gaussianFit(h_timing_thresh, fit_params_map[voltage].mean_timing_thresh, fit_params_map[voltage].sigma_timing_thresh);

        dir_amp->cd();       h_amplitude->Write();    delete h_amplitude;
        dir_rt->cd();        h_rise_time->Write();    delete h_rise_time;
        dir_fwhm->cd();      h_fwhm->Write();         delete h_fwhm;
        dir_charge->cd();    h_charge->Write();       delete h_charge;
        dir_timing_const->cd();  h_timing_const->Write();  delete h_timing_const;
        dir_timing_thresh->cd(); h_timing_thresh->Write(); delete h_timing_thresh;

        outputFile->cd();

        // 在循环中，处理完每个电压后，填充fit_tree
        voltage_var = voltage;
        mean_amp = fit_params_map[voltage].mean_amp;
        sigma_amp = fit_params_map[voltage].sigma_amp;
        mean_rt = fit_params_map[voltage].mean_rt;
        sigma_rt = fit_params_map[voltage].sigma_rt;
        mean_fwhm = fit_params_map[voltage].mean_fwhm;
        sigma_fwhm = fit_params_map[voltage].sigma_fwhm;
        mean_charge = fit_params_map[voltage].mean_charge;
        sigma_charge = fit_params_map[voltage].sigma_charge;
        mean_timing_CFD = fit_params_map[voltage].mean_timing_CFD;
        sigma_timing_CFD = fit_params_map[voltage].sigma_timing_CFD;
        mean_timing_thresh = fit_params_map[voltage].mean_timing_thresh;
        sigma_timing_thresh = fit_params_map[voltage].sigma_timing_thresh;

        fit_tree->Fill();



        cout << "电压 " << voltage << "V 拟合结果:" << endl;
        cout << "  幅值: " << fit_params_map[voltage].mean_amp << " ± " << fit_params_map[voltage].sigma_amp << " mV" << endl;
        cout << "  上升时间: " << fit_params_map[voltage].mean_rt << " ± " << fit_params_map[voltage].sigma_rt << " ns" << endl;
        cout << "  半高宽度: " << fit_params_map[voltage].mean_fwhm << " ± " << fit_params_map[voltage].sigma_fwhm << " ns" << endl;
        cout << "  电荷量: " << fit_params_map[voltage].mean_charge << " ± " << fit_params_map[voltage].sigma_charge << " e" << endl;
        cout << "  恒比定时差: " << fit_params_map[voltage].mean_timing_CFD << " ± " << fit_params_map[voltage].sigma_timing_CFD << " ns" << endl;
        cout << "  过阈定时差: " << fit_params_map[voltage].mean_timing_thresh << " ± " << fit_params_map[voltage].sigma_timing_thresh << " ns" << endl;
    }

    // 写入fit_tree并切换回根目录
    fit_tree->Write();
    outputFile->cd();

    // 创建二维图
    create2DPlotsOptimized(voltages, fit_params_map, amplitudes, rise_times, fwhms, charges,
                          timing_diffs_const, outputFile);

    // ======== TA修正和TQ修正 ========
    map<int, double> cfd_sigma_before, cfd_sigmaErr_before;
    map<int, double> cfd_sigma_afterTA, cfd_sigmaErr_afterTA;
    map<int, double> cfd_sigma_afterTQ, cfd_sigmaErr_afterTQ;
    map<int, double> thresh_sigma_before, thresh_sigmaErr_before;
    map<int, double> thresh_sigma_afterTA, thresh_sigmaErr_afterTA;
    map<int, double> thresh_sigma_afterTQ, thresh_sigmaErr_afterTQ;

    //
    performTA_TQ_Correction(voltages, fit_params_map, amplitudes, charges, rise_times, fwhms,
                            timing_diffs_const, timing_diffs_thresh, outputFile,
                            cfd_sigma_before, cfd_sigmaErr_before,
                            cfd_sigma_afterTA, cfd_sigmaErr_afterTA,
                            cfd_sigma_afterTQ, cfd_sigmaErr_afterTQ,
                            thresh_sigma_before, thresh_sigmaErr_before,
                            thresh_sigma_afterTA, thresh_sigmaErr_afterTA,
                            thresh_sigma_afterTQ, thresh_sigmaErr_afterTQ);

    // 输出汇总表
    cout << "\n\n========== 时间分辨汇总 ==========" << endl;
    cout << "电压(V)\tCFD_before(ns)\tCFD_TA(ns)\tCFD_TQ(ns)\tThresh_before(ns)\tThresh_TA(ns)\tThresh_TQ(ns)" << endl;
    for (int v : voltages) {
        cout << v << "\t"
             << (cfd_sigma_before.count(v) ? cfd_sigma_before[v] : 0) << " +/- "
             << (cfd_sigmaErr_before.count(v) ? cfd_sigmaErr_before[v] : 0) << "\t"
             << (cfd_sigma_afterTA.count(v) ? cfd_sigma_afterTA[v] : 0) << " +/- "
             << (cfd_sigmaErr_afterTA.count(v) ? cfd_sigmaErr_afterTA[v] : 0) << "\t"
             << (cfd_sigma_afterTQ.count(v) ? cfd_sigma_afterTQ[v] : 0) << " +/- "
             << (cfd_sigmaErr_afterTQ.count(v) ? cfd_sigmaErr_afterTQ[v] : 0) << "\t"
             << (thresh_sigma_before.count(v) ? thresh_sigma_before[v] : 0) << " +/- "
             << (thresh_sigmaErr_before.count(v) ? thresh_sigmaErr_before[v] : 0) << "\t"
             << (thresh_sigma_afterTA.count(v) ? thresh_sigma_afterTA[v] : 0) << " +/- "
             << (thresh_sigmaErr_afterTA.count(v) ? thresh_sigmaErr_afterTA[v] : 0) << "\t"
             << (thresh_sigma_afterTQ.count(v) ? thresh_sigma_afterTQ[v] : 0) << " +/- "
             << (thresh_sigmaErr_afterTQ.count(v) ? thresh_sigmaErr_afterTQ[v] : 0)
             << endl;

        // 保存修正后的时间分辨率到TAC_tree
        outputFile->cd();
        sigma_cfd_before = cfd_sigma_before.count(v) ? cfd_sigma_before[v] : 0;
        sigmaErr_cfd_before = cfd_sigmaErr_before.count(v) ? cfd_sigmaErr_before[v] : 0;
        sigma_cfd_afterTA = cfd_sigma_afterTA.count(v) ? cfd_sigma_afterTA[v] : 0;
        sigmaErr_cfd_afterTA = cfd_sigmaErr_afterTA.count(v) ? cfd_sigmaErr_afterTA[v] : 0;
        sigma_cfd_afterTQ = cfd_sigma_afterTQ.count(v) ? cfd_sigma_afterTQ[v] : 0;
        sigmaErr_cfd_afterTQ = cfd_sigmaErr_afterTQ.count(v) ? cfd_sigmaErr_afterTQ[v] : 0;
        sigma_thresh_before = thresh_sigma_before.count(v) ? thresh_sigma_before[v] : 0;
        sigmaErr_thresh_before = thresh_sigmaErr_before.count(v) ? thresh_sigmaErr_before[v] : 0;
        sigma_thresh_afterTA = thresh_sigma_afterTA.count(v) ? thresh_sigma_afterTA[v] : 0;
        sigmaErr_thresh_afterTA = thresh_sigmaErr_afterTA.count(v) ? thresh_sigmaErr_afterTA[v] : 0;
        sigma_thresh_afterTQ = thresh_sigma_afterTQ.count(v) ? thresh_sigma_afterTQ[v] : 0;
        sigmaErr_thresh_afterTQ = thresh_sigmaErr_afterTQ.count(v) ? thresh_sigmaErr_afterTQ[v] : 0;
        TAC_tree->Fill();
    }
    TAC_tree->Write();
    

    outputFile->Close();
    cout << "\n分析完成！结果保存在 " << outputFileName << endl;
    cout << "总共处理文件数: " << total_files << endl;
}

// 主函数
void RPC_PMT_DataAnalysis() {
    auto start_time = chrono::high_resolution_clock::now();

    analyzeDataSerial();

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
    cout << "Total analysis time: " << duration / 60. << " minutes" << endl;
}
