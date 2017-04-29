//
//  De-multiplex lane-based FASTQ files to sample-based FASTQ files.
//
//  main.cpp
//  Read-FASTQ-Files
//
//  Created by Yuguang Xiong on 4/26/17.
//  Copyright © 2017 Yuguang Xiong. All rights reserved.
//

#include <memory>
#include <array>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <cstdlib>
#include <csv.h>

enum SeqFragLine {Identifier=1, Sequence, Option, Quality};
const std::size_t n_seq_frag_lines = SeqFragLine::Quality;
// 2^17 = 128K
const std::size_t dge_conv_sample_buf_size = 131072;
const std::size_t dge_conv_sample_seq_frags_buf_size = dge_conv_sample_buf_size * n_seq_frag_lines;
const std::size_t well_bardcode_len = 6;
std::map<const std::string, const std::string> dge_conv_wells_samples = {{"ACTATT", "A8"}, {"GTAAAT", "B11"}, {"GTTATT", "B12"}, {"TCAATT", "C9"}, {"TCTATA", "C10"}, {"TGATAA", "C11"}, {"TGTTAT", "C12"}, {"TTATGT", "D3"}, {"TTCTAA", "D4"}, {"TCGAAG", "E3"}, {"GGACCA", "E5"}, {"ACCGCG", "E6"}, {"CCGCGA", "F3"}, {"CCGGTG", "F4"}, {"CGGGAG", "F11"}, {"CGGTCC", "F12"}};

int main(int argc, const char * argv[])
{
    if(argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " [FASTQ R1 file] [FASTQ R2 file]" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::string dge_fastq_file_path_r1(argv[1]);
    std::string dge_fastq_file_path_r2(argv[2]);
    
    int exit_code = EXIT_SUCCESS;
    
    try
    {
        bool check_fastq_file_lines = false;
        // Get the total number of lines of FASTQ file.
        if(check_fastq_file_lines)
        {
            io::LineReader dge_fastq_file_r1(dge_fastq_file_path_r1);
            while(dge_fastq_file_r1.next_line());
            auto dge_fastq_lines_r1 = dge_fastq_file_r1.get_file_line();
            io::LineReader dge_fastq_file_r2(dge_fastq_file_path_r2);
            while(dge_fastq_file_r2.next_line());
            auto dge_fastq_lines_r2 = dge_fastq_file_r2.get_file_line();
            
            if(dge_fastq_lines_r1 != dge_fastq_lines_r2)
            {
                std::cerr << "R1 lines = " << dge_fastq_lines_r1 << " and R2 lines = " << dge_fastq_lines_r2 << std::endl;
                std::cerr << "R1 and R2 file lines are different!" << std::endl;
                return EXIT_FAILURE;
            }
        }
        
        // Display the well barcodes of fragment sequences.
        //{
        // Initialize well barcodes.
        std::set<const std::string> dge_conv_well_barcodes;
        for(const auto& dge_conv_well_sample : dge_conv_wells_samples)
        {
            dge_conv_well_barcodes.insert(dge_conv_well_sample.first);
        }
        
        // Open lane-based input FASTQ R1 and R2 files.
        io::LineReader dge_fastq_file_r1(dge_fastq_file_path_r1);
        io::LineReader dge_fastq_file_r2(dge_fastq_file_path_r2);
        
        // Generate sample-based output FASTQ R1 and R2 file names.
        std::map<std::string, std::string> dge_conv_sample_fastq_file_paths_r1;
        std::map<std::string, std::string> dge_conv_sample_fastq_file_paths_r2;
        for(const auto& dge_conv_well_sample : dge_conv_wells_samples)
        {
            // Generate the name of DGE sample FASTQ file.
            // R1
            const auto pos_1 = dge_fastq_file_path_r1.find_last_of("\\/");
            std::string dge_fastq_file_dir_r1, dge_fastq_file_name_r1;
            if(pos_1 != std::string::npos)
            {
                dge_fastq_file_dir_r1 = dge_fastq_file_path_r1.substr(0, pos_1+1);
                dge_fastq_file_name_r1 = dge_fastq_file_path_r1.substr(pos_1+1);
            }
            else dge_fastq_file_name_r1 = dge_fastq_file_path_r1;
            std::string dge_conv_sample_fastq_file_path_r1 = dge_conv_well_sample.second + '.' + dge_fastq_file_name_r1 + ".txt";
            if(!dge_fastq_file_dir_r1.empty()) dge_conv_sample_fastq_file_path_r1 = dge_fastq_file_dir_r1 + dge_conv_sample_fastq_file_path_r1;
            dge_conv_sample_fastq_file_paths_r1[dge_conv_well_sample.first] = dge_conv_sample_fastq_file_path_r1;
            // R2
            const auto pos_2 = dge_fastq_file_path_r2.find_last_of("\\/");
            std::string dge_fastq_file_dir_r2, dge_fastq_file_name_r2;
            if(pos_2 != std::string::npos)
            {
                dge_fastq_file_dir_r2 = dge_fastq_file_path_r2.substr(0, pos_2+1);
                dge_fastq_file_name_r2 = dge_fastq_file_path_r2.substr(pos_2+1);
            }
            else dge_fastq_file_name_r2 = dge_fastq_file_path_r2;
            std::string dge_conv_sample_fastq_file_path_r2 = dge_conv_well_sample.second + '.' + dge_fastq_file_name_r2 + ".txt";
            if(!dge_fastq_file_dir_r2.empty()) dge_conv_sample_fastq_file_path_r2 = dge_fastq_file_dir_r2 + dge_conv_sample_fastq_file_path_r2;
            dge_conv_sample_fastq_file_paths_r2[dge_conv_well_sample.first] = dge_conv_sample_fastq_file_path_r2;
        }
        // Initialize sample-based output FASTQ R1 and R2 files.
        std::map<std::string, std::unique_ptr<std::ofstream>> dge_conv_sample_fastq_files_r1;
        for(const auto& dge_conv_sample_fastq_file_path_r1 : dge_conv_sample_fastq_file_paths_r1)
        {
            dge_conv_sample_fastq_files_r1[dge_conv_sample_fastq_file_path_r1.first] = std::unique_ptr<std::ofstream>(new std::ofstream(dge_conv_sample_fastq_file_path_r1.second));
        }
        std::map<std::string, std::unique_ptr<std::ofstream>> dge_conv_sample_fastq_files_r2;
        for(const auto& dge_conv_sample_fastq_file_path_r2 : dge_conv_sample_fastq_file_paths_r2)
        {
            dge_conv_sample_fastq_files_r2[dge_conv_sample_fastq_file_path_r2.first] = std::unique_ptr<std::ofstream>(new std::ofstream(dge_conv_sample_fastq_file_path_r2.second));
        }
        
        // Initialize sample-based FASTQ R1 and R2 sequence fragment buffers.
        std::map<std::string, std::vector<std::string>> dge_conv_sample_seq_frags_r1, dge_conv_sample_seq_frags_r2;
        for(const auto& dge_conv_well_sample : dge_conv_wells_samples)
        {
            dge_conv_sample_seq_frags_r1[dge_conv_well_sample.first] = std::vector<std::string>();
            dge_conv_sample_seq_frags_r2[dge_conv_well_sample.first] = std::vector<std::string>();
        }
        
        // Read R1 and R2 FASTQ files.
        std::array<std::string, n_seq_frag_lines> seq_frag_r1;
        std::array<std::string, n_seq_frag_lines> seq_frag_r2;
        bool found = false;
        std::string matched_barcode;
        std::size_t counts = 0;
        char* line_r1 = dge_fastq_file_r1.next_line();
        char* line_r2 = dge_fastq_file_r2.next_line();
        while(line_r1!=0 && line_r2!=0)
        {
            // Store each line of a fragment block.
            seq_frag_r1[counts] = line_r1;
            seq_frag_r2[counts] = line_r2;
            
            counts++;
            
            // At the 2nd line of a R1 fragment block that includes a well barcode.
            if(counts == SeqFragLine::Sequence)
            {
                const std::size_t line_len = std::strlen(line_r1);
                if(line_len >= well_bardcode_len)
                {
                    std::string barcode = ::strndup(line_r1, well_bardcode_len);
                    auto search = dge_conv_well_barcodes.find(barcode);
                    if(search != dge_conv_well_barcodes.end())
                    {
                        found = true;
                        matched_barcode = barcode;
                    }
                }
                else std::cerr << "Fragment sequence is too short!" << std::endl;
            }
            
            // At the end of a fragment block.
            if(counts == n_seq_frag_lines)
            {
                // Matched sequence fragment is found.
                if(found)
                {
                    std::vector<std::string>& dge_conv_sample_seq_frags_buf_r1 = dge_conv_sample_seq_frags_r1[matched_barcode];
                    std::vector<std::string>& dge_conv_sample_seq_frags_buf_r2 = dge_conv_sample_seq_frags_r2[matched_barcode];
                    
                    // Record matched fragment block.
                    for(const auto& seq_frag_line : seq_frag_r1) dge_conv_sample_seq_frags_buf_r1.push_back(seq_frag_line);
                    for(const auto& seq_frag_line : seq_frag_r2) dge_conv_sample_seq_frags_buf_r2.push_back(seq_frag_line);
                    
                    // Output when bufs are full.
                    if(dge_conv_sample_seq_frags_buf_r1.size() >= dge_conv_sample_seq_frags_buf_size)
                    {
                        for(const auto& seq_frag_line : dge_conv_sample_seq_frags_buf_r1) *dge_conv_sample_fastq_files_r1[matched_barcode] << seq_frag_line << std::endl;
                        dge_conv_sample_seq_frags_buf_r1.clear();
                    }
                    if(dge_conv_sample_seq_frags_buf_r2.size() >= dge_conv_sample_seq_frags_buf_size)
                    {
                        for(const auto& seq_frag_line : dge_conv_sample_seq_frags_buf_r2) *dge_conv_sample_fastq_files_r2[matched_barcode] << seq_frag_line << std::endl;
                        dge_conv_sample_seq_frags_buf_r2.clear();
                    }
                    // Reset found.
                    found = false;
                }
                counts = 0;
            }
            
            // Read next line.
            line_r1 = dge_fastq_file_r1.next_line();
            line_r2 = dge_fastq_file_r2.next_line();
        }
        
        // Output remaining bufs.
        for(const auto& dge_conv_sample_seq_frag_r1 : dge_conv_sample_seq_frags_r1)
        {
            for(const auto& seq_frag_line : dge_conv_sample_seq_frag_r1.second) *dge_conv_sample_fastq_files_r1[dge_conv_sample_seq_frag_r1.first] << seq_frag_line << std::endl;
        }
        for(const auto& dge_conv_sample_seq_frag_r2 : dge_conv_sample_seq_frags_r2)
        {
            for(const auto& seq_frag_line : dge_conv_sample_seq_frag_r2.second) *dge_conv_sample_fastq_files_r2[dge_conv_sample_seq_frag_r2.first] << seq_frag_line << std::endl;
        }
        //}
    }
    catch(std::exception& e)
    {
        exit_code = EXIT_FAILURE;
        std::cerr << "Error: " << e.what() << std::endl;
    }
    
    return exit_code;
}
