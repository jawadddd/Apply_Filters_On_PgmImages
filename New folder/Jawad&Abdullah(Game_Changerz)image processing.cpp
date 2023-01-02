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
        int row;
        row=Two.Rows;
        if(Rows>Two.Rows)
        {
            row=Rows;
        }


        for(int i=0;i<Two.Rows && i< MaxRows ;i++)
        {
            for(int j=Cols,k=0;j<Cols+Two.Cols && j<MaxCols;j++,k++)
            {
                Image[i][j]=Two.Image[i][k];
            }
        }
        if(Rows>Two.Rows)
        {
            int Left,Top,Bottom,Right;
             Left=Cols+1;
			 Top=Two.Rows+1;
			 Bottom=Rows;
			 Right=MaxCols;
		    Fill(Left,Top,Right,Bottom,fillValue);
	   	}
		if(Two.Rows>Rows)
		{
		    int Left,Top,Bottom,Right;
            Top=Rows+1;
			Bottom=Two.Rows;
			Left=0;
			Right=Cols;
	    	Fill(Left,Top,Right,Bottom,fillValue);
		}

		Rows=row;
        Cols=Cols+Two.Cols;

        return;
    }

    void combineTopToBottom(grayImage& Two, int fillValue = 0) {
        int col;
        col=Two.Cols;
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
            int Top,Bottom,Left,Right;
            Top=Two.Rows+1;
            Bottom=Rows;
            Left=Two.Cols;
            Right=Cols;
		    Fill(Left,Top,Right,Bottom,fillValue);
	   	}
		if(Cols<Two.Cols)
		{
			int Left,Top,Bottom,Right;
			 Top=Rows;
			 Bottom=Two.Rows;
			 Left=Cols;
			 Right=Two.Cols;
	    	Fill(Left,Top,Right,Bottom,fillValue);
		}
        Cols=col;
        Rows=Rows+Two.Rows;
        return;
    }

    void FadeIn(grayImage& Second,int frames, int seconds,string BaseFileName){

        grayImage G;
        double SizeX = 1.0/(frames*seconds);

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
        for(double Alpha = 1;Alpha>=0;Alpha-=SizeX)
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
                if(0>Alpha-SizeX && Alpha-SizeX > -SizeX)
                     Alpha = 0 ;
    }
    }

    int load(string File_Name){
        int columns;
        int rows;
        int MaxValue;
        string P2;
        string Comment;
        ifstream file(File_Name.c_str());
        getline(file,P2);
        getline(file,Comment);
        file>>columns>>rows;
        file>>Maximum;
        setRows(rows);
        setCols(columns);
        for(int i=0;i<Rows;i++)
            for(int j=0;j<Cols;j++)
                file>>Image[i][j];
        file.close();
        return 0;

    }

    int Save(string File_Name){

        ofstream file2(File_Name.c_str());

        file2<<"P2";
        file2<<endl<<"# Created by GAME CHANGERS "<<endl;
        file2<<Cols<<" "<<Rows<<endl;
        file2<<Maximum<<endl;

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

    void Rotate(grayImage& Rotated, double angle = 90, int aboutx = 0, int abouty = 0 )
    {
            int aa;

            int ab;

            angle=(3.14*angle)/180.0;

            for(int i=0;i<Rows;i++)
            {
                for(int j=0;j<Cols;j++)
                {
                    aa=((i-aboutx)*cos(angle)+(j-abouty)*sin(angle))+aboutx;
                    ab=((j-abouty)*cos(angle)-(i-aboutx)*sin(angle))+abouty;
                    Rotated.Image[i][j]=Image[aa][ab];
                }
            }

            Rotated.Rows=Rows;
            Rotated.Cols=Cols;
            Rotated.Maximum=Maximum;

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

    void Negative()
    {
        for(int i=0;i<Rows;i++)
            for(int j=0;j<Cols;j++)
                Image[i][j]=Maximum-Image[i][j];
    }

    void changeBrightness(int amount)
    {

        for(int i=0;i<Rows;i++)
        {
            for(int j=0;j<Cols;j++)
            {
                Image[i][j]=Image[i][j]+amount;



                if(Image[i][j]>Maximum)
                    Image[i][j]=Maximum;

                if(Image[i][j]<0)
                    Image[i][j]=0;


            }
        }
    }

    void Quantize(int Levels)
    {
        int p;
        p=Maximum/Levels;
        int A[Levels];
        cout<<endl;
        for(int a=0;a<Levels;a++)
        {
            cout<<"please Enter Number from "<<a*p<<" to "<<p*(a+1)<<" : "<<endl;
            cin>>A[a];
        }
        for(int a=0;a<Rows;a++)
            for(int b=0;b<Cols;b++)
                for(int c=0;c<Levels;c++)
                    if(Image[a][b]>c*p && Image[a][b]<(c+1)*p)
                        Image[a][b]=A[c];
        return;
    }

    void medianFilter(grayImage& Result, int filterSize = 3)
    {
        int A[filterSize*filterSize];
        for(int i=1;i<Rows-1;i++)
        {
            for(int j=1;j<Cols-1;j++)
            {
                int a=0;
                for(int k=i-(filterSize/2);k<=i+(filterSize/2);k++)
                {
                    for(int l=j-1;l<=j+1;l++)
                    {
                       A[a]=Image[k][l];
                       a++;
                    }
                }
                int b;
                for(int i=0;i<filterSize*filterSize;i++)
                {
                    for(int j=i+1;j<filterSize*filterSize;j++)
                    {
                        if(A[i]<A[j])
                        {
                            b=A[i];
                            A[i]=A[j];
                            A[j]=b;
                        }
                    }
                }
                Result.Image[i][j]=A[(filterSize*filterSize)/2];
            }
        }
        Result.Cols=Cols;
        Result.Rows=Rows;
    }
    void meanFilter(grayImage& Result, int filterSize = 3)
    {

        int added;
        for(int i=1;i<Rows-1;i++)
        {
            for(int j=1;j<Cols-1;j++)
            {
                added=0;
                for(int k=i-(filterSize/2);k<=i+(filterSize/2);k++)
                {
                    for(int l=j-1;l<=j+1;l++)
                    {
                       added=added+Image[k][l];
                    }
                }
                int avg=added/9;
                Result.Image[i][j]=avg;
            }
        }
        Result.Cols=Cols;
        Result.Rows=Rows;
    }
    void Resize(grayImage& Result, int NewRows, int NewColumns){
        double x,y;
        double m,n;
        x=NewRows;
        y=NewColumns;
        m=Rows/x;
        n=Cols/y;
        int k,l;
        int a;
        a=0;
        for(double i=0;i<Rows;i+=m)
        {
            int b;
            b=0;
            for(double j=0;j<Cols;j+=n)
            {
                k=i+0.9;
                l=j+0.9;

                Result.Image[a][b]=Image[k][l];
                b++;
            }
            a++;
        }
        Cols=NewColumns;
        Rows=NewRows;

    }
    void Resize(grayImage& Result,double Ratio){
        double Ro;
        Ro=Rows/Ratio;
        double Co;
        Co=Cols/Ratio;
        Resize(Result,Ro,Co);
    }
    void Transform(grayImage& RESULTANT,int Matrix[3][3])
    {
        double I;
        double J;
        double Z;
        for(int i=1;i<Rows;i++)
        {

            for(int j=1;j<Cols;j++)
            {
                RESULTANT.Image[((Matrix[0][0]*i)+(Matrix[0][1]*j)+(Matrix[0][2]))/((Matrix[2][0]*i)+(Matrix[2][1]*j)+(Matrix[2][2]))]
                    [((Matrix[1][0]*i)+(Matrix[1][1]*j)+(Matrix[1][2]))/((Matrix[2][0]*i)+(Matrix[2][1]*j)+(Matrix[2][2]))]=Image[i][j];
            }

        }

    }
    void Filter(double Mask[3][3])
    {
        int a;
        int sum;
        int b;
        for(int l=1;l<Rows-1;l++)
        {
            for(int k=1;k<Cols-1;k++)
            {
                 a=0;
                 sum=0;
                for(int i=l-1;i<=l+1;i++)
                {
                    b=0;
                    for(int j=k-1;j<=k+1;j++)
                    {
                        sum=sum+(Mask[a][b]*Image[i][j]);
                        b++;
                    }
                    a++;
                }



                Image[l][k]=sum;
            }
        }
    }
    void DerivativeImage(){
        int Mask[3][3]={{-1,0,1},{-1,0,1},{-1,0,1}};
        int Maskx[3][3]={{-1,-1,-1},{0,0,0},{1,1,1}};
        int a;
        int b;
        int sum2;
        int sq;
        for(int l=1;l<Rows-1;l++)
        {
            for(int k=1;k<Cols-1;k++)
            {
                a=0;
                int sum1=0;
                for(int i=l-1;i<=l+1;i++)
                {
                    int t=0;
                    for(int j=k-1;j<=k+1;j++)
                    {
                        sum1=sum1+(Mask[a][t]*Image[i][j]);
                        t++;
                    }
                    a++;
                }
                b=0;
                sum2=0;
                for(int i=l-1;i<=l+1;i++)
                {
                    int t=0;
                    for(int j=k-1;j<=k+1;j++)
                    {
                        sum2=sum2+(Maskx[b][t]*Image[i][j]);
                        t++;
                    }
                    b++;
                }
                sq=sqrt((sum1/9)*(sum1/9)+(sum2/9)*(sum2/9));
                Image[l][k]=sq;
            }
        }

    }

    void Crop(int top,int bottom,int left,int right){

        for(int i=top,xi=0;i<bottom;i++,xi++)
            for(int j=left,xj=0;j<right;j++,xj++)
                Image[xi][xj]=Image[i][j];

        Cols=right-left;
        Rows=bottom-top;
    }

private:

    void FilpHorizontal()
    {
        int temporary[MaxRows][MaxCols];
        for(int i=Rows-1,xi=0;i>=0 && xi<Rows ;i--,xi++)
            for(int j=0;j<Cols;j++)
                temporary[xi][j]=Image[i][j];
        for(int xi=0; xi<Rows ;xi++)
            for(int j=0;j<Cols;j++)
                Image[xi][j]=temporary[xi][j];
    }

    void FlipVerttical()
    {
        int temporary[MaxRows][MaxCols];
        for(int i=0;i<Rows;i++)
            for(int j=0,xi=Cols;j<Cols;xi--,j++)
                temporary[i][j]=Image[i][xi];
        for(int i=0;i<Rows;i++)
            for(int j=0;j<Cols;j++)
                Image[i][j]=temporary[i][j];
    }

    void Fill(int L, int T, int R, int B, int FillValue)
    {
        for(int i = T; i<= B; i++)
            for(int j = L; j <= R; j++)
                Image[i][j] = FillValue;
    }

    unsigned short Image[MaxRows][MaxCols];
    int Rows;
    int Cols;
    int Maximum;
    bool loaded;
};

int main()
{
    grayImage GM, GM2;
    cout<<"         (MENU)             "<<endl;
    cout<<"1.Combine Top to Bottom "<<endl<<"2.Combine Side by Side "<<endl<<"3.Flip "<<endl<<"4.Change Brightness "<<endl<<"5.Negative "
    <<endl<<"6.Quantize of an Image "<<endl<<"7.Median Filter "<<endl<<"8.Mean Filter "<<endl<<"9.Filter(MASK) "<<endl<<"10.Resize by NEW ROWS and COLUMNS "<<endl<<"11.Resize by Ratio"
    <<endl<<"12.Rotation "<<endl<<"13.Transform of an Image "<<endl<<"14.Derivative of an Image "<<endl<<"15.Crop "<<endl<<"16.FADEIN";

    cout<<endl<<"Choose the number on basis of the menu plz!! ";
    cout<<endl;
    int xyz;
    cin>>xyz;
    string FIL1;
    string FIL2;
    string FIL3;
    cout<<"Enter the FILE NAME for processing : ";
    cout<<endl;
    cin>>FIL1;
    GM.load(FIL1);
    if(xyz==1)
    {

        cout<<"Enter the 2nd FILE NAME for Combining two Images in form of top by bottom "<<endl;;
        cin>>FIL2;
        GM2.load(FIL2);
        GM.combineTopToBottom(GM2, 150);
    }
    if(xyz==2)
    {
        cout<<"Enter the 2nd FILE NAME for Combining two Images in form of side by side "<<endl;
        cin>>FIL2;
        GM2.load(FIL2);
        GM.combineSideBySide(GM2, 150);
    }
    if(xyz==3)
    {
        int Flip;
        cout<<"1.Vertical Flip "<<endl<<"2.Horizontal Flip "<<endl;
        cout<<"choose 1 or 2: "<<endl;
        cin>>Flip;
        if(Flip==1)
        {
            GM.Flip(1);
        }

        if(Flip==2)
        {
            GM.Flip();
        }
    }
    if(xyz==4)
    {
        cout<<endl<<"HOW DO YOU WANT TO CHANGE THE BRIGHTNESS OF YOUR ENTERED IMAGE!!ENTER ANY AMOUNT:"<<endl;
        int amount;
        cin>>amount;
        GM.changeBrightness(amount);
    }
    if(xyz==5)
    {

        cout<<endl;
        GM.Negative();

    }
    if(xyz==6)
    {
        int Levels;
        cout<<"Enter the Number of Levels: ";
        cin>>Levels;
        GM.Quantize(Levels);
    }
    if(xyz==7)
    {

        GM.medianFilter(GM,3);
    }
    if(xyz==8)
    {

        GM.meanFilter(GM,3);
    }
    if(xyz==9)
    {
        double Mask[3][3]={{0.15,0.15,0.15},{0.15,0.15,0.15},{0.15,0.15,0.15}};
        GM.Filter(Mask);
    }
    if(xyz==10)
    {
        int modifiedRows;
        int modifiedColumns;
        cout<<"Enter the new Rows number: "<<endl;
        cin>>modifiedRows;
        cout<<"Enter the new Columns number: "<<endl;
        cin>>modifiedColumns;
        GM.Resize(GM,modifiedRows,modifiedColumns);
    }
    if(xyz==11)
    {
        double Ratio;
        cout<<"Enter Ratio plz!! "<<endl;
        cin>>Ratio;
        GM.Resize(GM,Ratio);
    }
    if(xyz==12)
    {
        int Angle;
        int x;
        int y;
        cout<<"Enter Angle plz!! "<<endl;
        cin>>Angle;
        cout<<"Enter Angle About x-axis: "<<endl;
        cin>>x;
        cout<<"Enter Angle About y-axis: "<<endl;
        cin>>y;

        GM.Rotate(GM,Angle,x,y);
    }
    if(xyz==13)
    {
        int MASK[3][3]={{1,0,0},{1,-1,1},{0,0,1}};
        GM.Transform(GM,MASK);
    }
    if(xyz==14)
    {
        cout<<endl;
        GM.DerivativeImage();
    }
    if(xyz==15)
    {
        int top;
        int bottom;
        int left;
        int right;
        cout<<"Enter Top plz!!! "<<endl;
        cin>>top;
        cout<<"Enter Bottom plz!! "<<endl;
        cin>>bottom;
        cout<<"Enter Left plz!! "<<endl;
        cin>>left;
        cout<<"Enter Right plz!! "<<endl;
        cin>>right;
        GM.Crop(top,bottom,left,right);
    }
    if(xyz==16)
    {
        int Frames,Seconds;
        cout<<"Enter Frames plz!! "<<endl;
        cin>>Frames;
        cout<<"Enter Seconds plz!! "<<endl;
        cin>>Seconds;
        cout<<"Enter the name of second file 4 Fade_in!!! ";
        cout<<endl;
        cin>>FIL2;
        GM2.load(FIL2);
        GM.FadeIn(GM2, Frames, Seconds, "Results\\Second\\");
        return 0;
    }
    cout<<endl;
    cout<<"Enter new filename for saving Data: ";
    cout<<endl;
    cin>>FIL3;
    GM.Save(FIL3);
    return 0;
}
