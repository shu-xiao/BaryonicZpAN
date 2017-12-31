#include <iostream>
#include <fstream>

using namespace std;
int* nEventList (string fileName) {
    int* data;
    const int size = 6;
    ifstream inputFile(fileName.data());
    data = new int[size];  
    for(int i=0;i<size;i++)
    {
        inputFile >> data[i];
    }
    return data;
}
void effiCutPlot() {
    int zpMass[8] = {600,800,1000,1200,1400,1700,2000,2500};
    vector <int*> eventList;
    eventList.reserve(8);
    for (int n=0;n<8;n++) {
        string fname = Form("effi_Zpmass%d_A0mass300_4jets.txt",zpMass[n]); 
        int* nEvent = nEventList(fname);
        /*
        for (int i=0;i<6;i++) {
            cout << *nEvent << endl;
            nEvent++;
        }
        */
        eventList.push_back(nEvent);
    }
    for (int i=0;i<8;i++) {
        for (int j=0;j<6;j++) {
            cout << *eventList.at(i) << "\t";
            eventList.at(i)++;
        }
        cout << endl;
    }

}
