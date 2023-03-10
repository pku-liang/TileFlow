#pragma once 

#define MCARO_WRAPPER(code) do {code} while(0) 

#define TILEFLOW_ERROR(msg) do{std::cerr << "[ERROR]: " << msg << std::endl;exit(1);} while(0)

#define TILEFLOW_ASSERT(cond, msg) do{if(!(cond)) {std::cerr << "[ASSERT ERROR]: " << msg << std::endl; exit(1);} }while(0)

#define TILEFLOW_WARNING(msg) do{std::cerr << "[WARNING]: " << msg << std::endl;}while(0)

const int MaxTensors = 32;