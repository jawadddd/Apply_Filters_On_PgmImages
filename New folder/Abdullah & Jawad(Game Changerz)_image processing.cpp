#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<cstdlib>

using namespace std;
bool Saved = false;

#define MaxRows 500
#define MaxCols 500

struct grayImage{

    grayImage(){
        Rows = Cols = 0;
        loaded = false;
        for(int r = 0; r< MaxRows; r++)
            for(int c = 0; c< MaxCols; c++)
                Image[r][c] = 0;
    }

    unsigned short setPixel(unsigned short value, int r, int c) {
        if (r >= Rows || c >= Cols || r < 0 || c < 0 || value < 0 || value > Maximum)
            return -1;

        Image[r][c] = value;
        return value;
    }

    int getPixel(int r, int c){
        if (r >= Rows || c >= Cols || r < 0 || c < 0)
            return -1;
        return Image[r][c];
    }

    int setRows(int rows){
        if (rows < 1 || rows > MaxRows)
            return -1;
        Rows = rows;
        return Rows;
    }

    int getRows(){
        return Rows;
    }

    int setCols(int cols){
        if (cols < 1 || cols > MaxCols)
            return -1;
        Cols = cols;
        return Cols;
    }

    int getCols() {
        return Cols;
    }

    void combineSideBySide(grayImage& Two, int fillValue = 0){
        int row=Two.Rows;
        if(Rows>Two.Rows)
            row=Rows;
        for(int i=0;i<Two.Rows && i< MaxRows ;i++)
        {
            for(int j=Cols,k=0;j<Cols+Two.Cols && j<MaxCols;j++,k++)
            {
                Image[i][j]=Two.Image[i][k];
            }
        }
        if(Rows>Two.Rows){
			int Top=Two.Rows+1;
			int Bottom=Rows;
			int Left=Cols+1;
			int Right=MaxCols;
		    Fill(Left,Top,Right,Bottom,fillValue);
	   	}
		if(Two.Rows>Rows)
		{
			int Top=Rows+1;
			int Bottom=Two.Rows;
			int Left=0;
			int Right=Cols;
	    	Fill(Left,Top,Right,Bottom,fillValue);
		}
        Cols=Cols+Two.Cols;
        Rows=row;
        return;
    }

    void combineTopToBottom(grayImage& Two, int fillValue = 0) {
        int col=Two.Cols;
        if(Cols>Two.Cols)
            col=Cols;
        for(int i=Rows,k=0;k<Two.Rows && i<MaxRows;i++,k++)
        {
            for(int j=0;j<Two.Cols && j<MaxCols;j++)
            {
                Image[i][j]=Two.Image[k][j];
            }
        }
        if(Cols>Two.Cols)
        {
			int Top=Two.Rows+1;
			int Bottom=Rows;
			int Left=Two.Cols;
			int Right=Cols;
		    Fill(Left,Top,Right,Bottom,fillValue);
	   	}
		if(Cols<Two.Cols)
		{
			int Top=Rows;
			int Bottom=Two.Rows;
			int Left=Cols;
			int Right=Two.Cols;
	    	Fill(Left,Top,Right,Bottom,fillValue);
		}
        Cols=col;
        Rows=Rows+Two.Rows;
        return;
    }

    void FadeIn(grayImage& Second,int frames, int seconds,string BaseFileName){

        grayImage G;
        double StepSize = 1.0/(frames*seconds);

        int R = Rows;
        if(R < Second.Rows)
            R = Second.Rows;
        int C = Cols;
        if(C < Second.Cols)
            C = Second.Cols;
        G.Rows = R;
        G.Cols = C;
        G.loaded = true;
        G.Maximum = Maximum;
        if(Maximum < Second.Maximum)
            G.Maximum = Second.Maximum;

        int Counter = 0;
        for(double Alpha = 1;Alpha>=0;Alpha-=StepSize)
        {
            for(int i = 0; i < R; i++)
            {
                for(int j = 0;j < C; j++)
                {
                    G.Image[i][j] = Image[i][j]*Alpha + (1-Alpha)*Second.Image[i][j];
                }
            }
            char N[10];
            itoa(Counter,N,10);
            G.Save(BaseFileName + N + ".pgm");
            Counter ++;
                if(0>Alpha-StepSize && Alpha-StepSize > -StepSize)
                     Alpha = 0 ;
    }
    }

    int load(string File_Name){
        int columns,rows,MaxValue;
        string P2,Comment;
        ifstream file(File_Name.c_str());
        getline(file,P2);
        getline(file,Comment);
        file>>columns>>rows>>Maximum;
        setRows(rows);
        setCols(columns);
        for(int i=0;i<Rows;i++)
        {

            for(int j=0;j<Cols;j++)
            {

                file>>Image[i][j];
            }
        }
        file.close();
        return 0;

    }

    int Save(string File_Name){

        ofstream file2(File_Name.c_str());
        file2<<"P2"<<endl<<"# Created by GAME CHANGERS "<<endl<<Cols<<" "<<Rows<<endl<<Maximum<<endl;
        for(int i=0;i<Rows;i++)
        {

            for(int j=0;j<Cols;j++)
            {

                file2<<Image[i][j]<<" ";
            }
            file2<<endl;
        }
        file2.close();
        return 0;

    }

    void Rotate(grayImage& Two, double angle = 90, int aboutx = 0, int abouty = 0 ){
            int m,k;
            angle=(3.14*angle)/180.0;
            for(int i=0;i<Rows;i++)
            {
                for(int j=0;j<Cols;j++)
                {
                    m=((i-aboutx)*cos(angle)+(j-abouty)*sin(angle))+aboutx;
                    k=((j-abouty)*cos(angle)-(i-aboutx)*sin(angle))+abouty;
                    Two.Image[i][j]=Image[m][k];
                }
            }
            Two.Rows=Rows;
            Two.Cols=Cols;
            Two.Maximum=Maximum;
    }

    void Flip(int HorizontalOrVertical = 0){
        if(HorizontalOrVertical == 1)
        {
            FlipVerttical();
        }
        else
        {
            FilpHorizontal();
        }
        return;

    }

    void Negative(){
        for(int i=0;i<Rows;i++)
        {
            for(int j=0;j<Cols;j++)
            {
                Image[i][j]=Maximum-Image[i][j];
            }
        }

    }

    void changeBrightness(int amount){
        for(int i=0;i<Rows;i++)
        {
            for(int j=0;j<Cols;j++)
            {
                Image[i][j]=Image[i][j]+amount;
                if(Image[i][j]>Maximum)
                {
                    Image[i][j]=Maximum;
                }
                if(Image[i][j]<0)
                {
                    Image[i][j]=0;
                }
            }
        }
    }

    void Quantize(int Levels){
        int p=Maximum/Levels;
        int A[Levels];
        cout<<endl;
        for(int i=0;i<Levels;i++)
        {
            cout<<"Enter the Number from "<<i*p<<" to "<<p*(i+1)<<" : ";
            cin>>A[i];
        }
        for(int i=0;i<Rows;i++)
        {
            for(int j=0;j<Cols;j++)
            {
                for(int k=0;k<Levels;k++)
                {
                    if(Image[i][j]>k*p && Image[i][j]<(k+1)*p)
                    {
                        Image[i][j]=A[k];
                    }
                }
            }
        }
        return;
    }

    void medianFilter(grayImage& Result, int filterSize = 3){
        int A[filterSize*filterSize];
        for(int i=1;i<Rows-1;i++)
        {
            for(int j=1;j<Cols-1;j++)
            {
                int b=0;
                for(int k=i-(filterSize/2);k<=i+(filterSize/2);k++)
                {
                    for(int p=j-1;p<=j+1;p++)
                    {
                       A[b]=Image[k][p];
                       b++;
                    }
                }
                int g;
                for(int i=0;i<filterSize*filterSize;i++)
                {
                    for(int j=i+1;j<filterSize*filterSize;j++)
                    {
                        if(A[i]<A[j])
                        {
                            g=A[i];
                            A[i]=A[j];
                            A[j]=g;
                        }
                    }
                }
                Result.Image[i][j]=A[(filterSize*filterSize)/2];
            }
        }
        Result.Rows=Rows;
        Result.Cols=Cols;
    }

    void meanFilter(grayImage& Result, int filterSize = 3){

        int sum,avg;
        for(int i=1;i<Rows-1;i++)
        {
            for(int j=1;j<Cols-1;j++)
            {
                sum=0;
                for(int k=i-(filterSize/2);k<=i+(filterSize/2);k++)
                {
                    for(int p=j-1;p<=j+1;p++)
                    {
                       sum=sum+Image[k][p];
                    }
                }
                avg=sum/9;
                Result.Image[i][j]=avg;
            }
        }
        Result.Rows=Rows;
        Result.Cols=Cols;
    }

    void Resize(grayImage& Result, int NewRows, int NewColumns){
        double x=NewRows,y=NewColumns,p,q;
        p=Rows/x;
        q=Cols/y;
        int k,t,e=0;
        for(double i=0;i<Rows;i+=p)
        {
            int f=0;
            for(double j=0;j<Cols;j+=q)
            {
                k=i+0.9;
                t=j+0.9;

                Result.Image[e][f]=Image[k][t];
                f++;
            }
            e++;
        }
        Rows=NewRows;
        Cols=NewColumns;
    }

    void Resize(grayImage& Result,double Ratio){
        double R=Rows/Ratio,C=Cols/Ratio;
        Resize(Result,R,C);
    }

    void Transform(grayImage& Two,int Mask[3][3]){
        double I,J,Z;
        for(int i=1;i<Rows;i++)
        {
            for(int j=1;j<Cols;j++)
            {
                Two.Image[((Mask[0][0]*i)+(Mask[0][1]*j)+(Mask[0][2]))/((Mask[2][0]*i)+(Mask[2][1]*j)+(Mask[2][2]))]
                    [((Mask[1][0]*i)+(Mask[1][1]*j)+(Mask[1][2]))/((Mask[2][0]*i)+(Mask[2][1]*j)+(Mask[2][2]))]=Image[i][j];
            }

        }
    }

    void Filter(double Mask[3][3]){
        for(int p=1;p<Rows-1;p++)
        {
            for(int k=1;k<Cols-1;k++)
            {
                int f=0;
                int sum=0;
                for(int i=p-1;i<=p+1;i++)
                {
                    int t=0;
                    for(int j=k-1;j<=k+1;j++)
                    {
                        sum=sum+(Mask[f][t]*Image[i][j]);
                        t++;
                    }
                    f++;
                }
                Image[p][k]=sum;
            }
        }
    }

    void DerivativeImage(){
        int Mask[3][3]={{-1,0,1},{-1,0,1},{-1,0,1}};
        int Mask2[3][3]={{-1,-1,-1},{0,0,0},{1,1,1}};
        int f,f1,sum2;
        for(int p=1;p<Rows-1;p++)
        {
            for(int k=1;k<Cols-1;k++)
            {
                f=0;
                int sum1=0;
                for(int i=p-1;i<=p+1;i++)
                {
                    int t=0;
                    for(int j=k-1;j<=k+1;j++)
                    {
                        sum1=sum1+(Mask[f][t]*Image[i][j]);
                        t++;
                    }
                    f++;
                }
                f1=0;
                sum2=0;
                for(int i=p-1;i<=p+1;i++)
                {
                    int t=0;
                    for(int j=k-1;j<=k+1;j++)
                    {
                        sum2=sum2+(Mask2[f1][t]*Image[i][j]);
                        t++;
                    }
                    f1++;
                }
                int sq=sqrt((sum1/9)*(sum1/9)+(sum2/9)*(sum2/9));
                Image[p][k]=sq;
            }
        }

    }

    void Crop(int top,int bottom,int left,int right){

        for(int i=top,a=0;i<bottom;i++,a++)
        {
            for(int j=left,b=0;j<right;j++,b++)
            {
                Image[a][b]=Image[i][j];
            }
        }
        Rows=bottom-top;
        Cols=right-left;

    }

private:

    void FilpHorizontal(){
        int temp[MaxRows][MaxCols];
        for(int i=Rows-1,d=0;i>=0 && d<Rows ;i--,d++)
        {
            for(int j=0;j<Cols;j++)
            {
                temp[d][j]=Image[i][j];
            }
        }
        for(int d=0; d<Rows ;d++)
        {
            for(int j=0;j<Cols;j++)
            {
                Image[d][j]=temp[d][j];
            }
        }
    }

    void FlipVerttical(){
        int temp[MaxRows][MaxCols];
        for(int i=0;i<Rows;i++)
        {
            for(int j=0,d=Cols;j<Cols;d--,j++)
            {
                temp[i][j]=Image[i][d];
            }
        }
        for(int i=0;i<Rows;i++)
        {
            for(int j=0;j<Cols;j++)
            {
                Image[i][j]=temp[i][j];
            }
        }
    }

    void Fill(int L, int T, int R, int B, int FillValue){
        for(int i = T; i<= B; i++)
            for(int j = L; j <= R; j++)
                Image[i][j] = FillValue;
    }

    unsigned short Image[MaxRows][MaxCols];
    int Rows, Cols, Maximum;
    bool loaded;
};

int main()
{
    grayImage GM, GM2;
    cout<<"                              ---------------MAIN MENU------------------ "<<endl;
    cout<<"1.Combine Top to Bottom "<<endl<<"2.Combine Side by Side "<<endl<<"3.Flip "<<endl<<"4.Change Brightness "<<endl<<"5.Negative "
    <<endl<<"6.Quantize of an Image "<<endl<<"7.Median Filter "<<endl<<"8.Mean Filter "<<endl<<"9.Filter(MASK) "<<endl<<"10.Resize by NEW ROWS and COLUMNS "<<endl<<"11.Resize by Ratio"
    <<endl<<"12.Rotation "<<endl<<"13.Transform of an Image "<<endl<<"14.Derivative of an Image "<<endl<<"15.Crop "<<endl<<"16.FADEIN";

    cout<<endl<<"Enter Your Choice: ";
    int n;
    cin>>n;
    string f1,f2,f3;
    cout<<"Enter the FILE NAME for processing : ";
    cin>>f1;
    GM.load(f1);
    if(n==1)
    {

        cout<<"Enter the 2nd FILE NAME for Combining two Images ";
        cin>>f2;
        GM2.load(f2);
        GM.combineTopToBottom(GM2, 150);
    }
    if(n==2)
    {
        cout<<"Enter the 2nd FILE NAME for Combining two Images ";
        cin>>f2;
        GM2.load(f2);
        GM.combineSideBySide(GM2, 150);
    }
    if(n==3)
    {
        int Flip;
        cout<<"1.Vertical Flip "<<endl<<"2.Horizontal Flip "<<endl;
        cout<<"Enter your Choice: ";
        cin>>Flip;
        if(Flip==1)
        {
            GM.Flip(1);
        }
        else
        {
            GM.Flip();
        }
    }
    if(n==4)
    {
        GM.changeBrightness(200);
    }
    if(n==5)
    {
        GM.Negative();
    }
    if(n==6)
    {
        int Levels;
        cout<<"Enter the Number of Levels: ";
        cin>>Levels;
        GM.Quantize(Levels);
    }
    if(n==7)
    {
        GM.medianFilter(GM,3);
    }
    if(n==8)
    {
        GM.meanFilter(GM,3);
    }
    if(n==9)
    {
        double Mask[3][3]={{0.15,0.15,0.15},{0.15,0.15,0.15},{0.15,0.15,0.15}};
        GM.Filter(Mask);
    }
    if(n==10)
    {
        int NewRows,NewColumns;
        cout<<"Enter NewRows: ";
        cin>>NewRows;
        cout<<"Enter NewColums: ";
        cin>>NewColumns;
        GM.Resize(GM,NewRows,NewColumns);
    }
    if(n==11)
    {
        double Ratio;
        cout<<"Enter the Ratio: ";
        cin>>Ratio;
        GM.Resize(GM,Ratio);
    }
    if(n==12)
    {
        int Angle,x,y;
        cout<<"Enter the Angle: ";
        cin>>Angle;
        cout<<"About x-axis: ";
        cin>>x;
        cout<<"About y-axis: ";
        cin>>y;

        GM.Rotate(GM,Angle,x,y);
    }
    if(n==13)
    {
        int MASK[3][3]={{1,0,0},{1,-1,1},{0,0,1}};
        GM.Transform(GM,MASK);
    }
    if(n==14)
    {
        GM.DerivativeImage();
    }
    if(n==15)
    {
        int top,bottom,left,right;
        cout<<"Enter the Top: ";
        cin>>top;
        cout<<"Enter the Bottom: ";
        cin>>bottom;
        cout<<"Enter the Left: ";
        cin>>left;
        cout<<"Enter the Right: ";
        cin>>right;
        GM.Crop(top,bottom,left,right);
    }
    if(n==16)
    {
        int Frames,Seconds;
        cout<<"Enter the Frames: ";
        cin>>Frames;
        cout<<"Enter the Seconds: ";
        cin>>Seconds;
        cout<<"Enter the 2nd FILE NAME for FadeIN ";
        cin>>f2;
        GM2.load(f2);
        GM.FadeIn(GM2, Frames, Seconds, "Results\\Second\\");
        return 0;
    }
    cout<<"Enter the FILE NAME for SAVING Data: ";
    cin>>f3;
    GM.Save(f3);
    return 0;
}
