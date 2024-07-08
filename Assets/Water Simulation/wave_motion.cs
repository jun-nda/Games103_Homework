using UnityEngine;
using System.Collections;

public class wave_motion : MonoBehaviour 
{
	int size 		= 100;
	float rate 		= 0.005f;
	float gamma		= 0.004f;
	float damping 	= 0.996f;
	float[,] 	old_h;
	float[,]	low_h;
	float[,]	vh;
	float[,]	b;

	bool [,]	cg_mask;
	float[,]	cg_p;
	float[,]	cg_r;
	float[,]	cg_Ap;
	bool 	tag = true;

	Vector3 	cube_v = Vector3.zero;
	Vector3 	cube_w = Vector3.zero;


	// Use this for initialization
	void Start () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.Clear ();

		Vector3[] X=new Vector3[size*size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			X[i*size+j].x=i*0.1f-size*0.05f;
			X[i*size+j].y=0;
			X[i*size+j].z=j*0.1f-size*0.05f;
		}

		int[] T = new int[(size - 1) * (size - 1) * 6];
		int index = 0;
		for (int i=0; i<size-1; i++)
		for (int j=0; j<size-1; j++)
		{
			T[index*6+0]=(i+0)*size+(j+0);
			T[index*6+1]=(i+0)*size+(j+1);
			T[index*6+2]=(i+1)*size+(j+1);
			T[index*6+3]=(i+0)*size+(j+0);
			T[index*6+4]=(i+1)*size+(j+1);
			T[index*6+5]=(i+1)*size+(j+0);
			index++;
		}
		mesh.vertices  = X;
		mesh.triangles = T;
		mesh.RecalculateNormals ();

		low_h 	= new float[size,size];
		old_h 	= new float[size,size];
		vh 	  	= new float[size,size];
		b 	  	= new float[size,size];

		cg_mask	= new bool [size,size];
		cg_p 	= new float[size,size];
		cg_r 	= new float[size,size];
		cg_Ap 	= new float[size,size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			low_h[i,j]=99999;
			old_h[i,j]=0;
			vh[i,j]=0;
		}
	}

	void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj)
	{
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			Ax[i,j]=0;
			if(i!=0)		Ax[i,j]-=x[i-1,j]-x[i,j];
			if(i!=size-1)	Ax[i,j]-=x[i+1,j]-x[i,j];
			if(j!=0)		Ax[i,j]-=x[i,j-1]-x[i,j];
			if(j!=size-1)	Ax[i,j]-=x[i,j+1]-x[i,j];
		}
	}

	float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj)
	{
		float ret=0;
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			ret+=x[i,j]*y[i,j];
		}
		return ret;
	}

	void Conjugate_Gradient(bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj)
	{
		//Solve the Laplacian problem by CG.
		A_Times(mask, x, cg_r, li, ui, lj, uj);

		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			cg_p[i,j]=cg_r[i,j]=b[i,j]-cg_r[i,j];
		}

		float rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);

		for(int k=0; k<128; k++)
		{
			if(rk_norm<1e-10f)	break;
			A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
			float alpha=rk_norm/Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				x[i,j]   +=alpha*cg_p[i,j];
				cg_r[i,j]-=alpha*cg_Ap[i,j];
			}

			float _rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);
			float beta=_rk_norm/rk_norm;
			rk_norm=_rk_norm;

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				cg_p[i,j]=cg_r[i,j]+beta*cg_p[i,j];
			}
		}

	}

	void Shallow_Wave(float[,] old_h, float[,] h, float [,] new_h)
	{		
		//Step 1:
		//TODO: Compute new_h based on the shallow wave model.
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				new_h[i, j] = h[i, j] + (h[i, j] - old_h[i, j]) * damping; // 先是时序上的，自己前一帧和前两帧的影响
				
				// 然后是空间上的，自己周围四个邻居对自己的影响
				if (i - 1 >= 0) new_h[i, j] += (h[i - 1, j] - h[i, j]) * rate;
				if (i + 1 < size) new_h[i, j] += (h[i + 1, j] - h[i, j]) * rate;
				if (j - 1 >= 0) new_h[i, j] += (h[i, j - 1] - h[i, j]) * rate;
				if (j + 1 < size) new_h[i, j] += (h[i, j + 1] - h[i, j]) * rate;
			}
		}
		
		//Step 2: Block->Water coupling
		//TODO: for block 1, calculate low_h.
		// 先找到方块
		GameObject cube = GameObject.Find("Block");
		Vector3 cubePos = cube.transform.position;
		Mesh cubeMesh = cube.GetComponent<MeshFilter>().mesh;
		
		// 参考答案里有个+-3不知道干啥的，去掉感觉效果更好点
		int lower_i = (int)((cubePos.x + 5.0f) * 10);
		int upper_i = (int)((cubePos.x + 5.0f) * 10);
		int lower_j = (int)((cubePos.z + 5.0f) * 10);
		int upper_j = (int)((cubePos.z + 5.0f) * 10);
		Bounds bounds = cubeMesh.bounds;
		
		for (int i = lower_i; i <= upper_i; i++){
			for (int j = lower_j; j <= upper_j; j++){
				if (i >= 0 && j >= 0 && i < size && j < size){
					Vector3 p = new Vector3(i * 0.1f - size * 0.05f, 0, j * 0.1f - size * 0.05f); // -5 到 +5， 咱也不知道写这么麻烦干嘛
					Vector3 q = new Vector3(i * 0.1f - size * 0.05f, 1, j * 0.1f - size * 0.05f);
					p = cube.transform.InverseTransformPoint(p); // 世界转局部
					q = cube.transform.InverseTransformPoint(q);

					// 在立方体自己的空间里，算一个相交的高度
					Ray ray = new Ray(p, q - p);
					bounds.IntersectRay(ray, out var dist);

					low_h[i, j] =  dist;//cube_p.y-0.5f;  没有双向影响得话，就是-0.5f
				}
			}
		}
		
		//TODO: then set up b and cg_mask for conjugate gradient.
		for (int i = 0; i < size; i++){
			for (int j = 0; j < size; j++)
			{
				if (low_h[i, j] > h[i, j])
				{
					b[i, j] = 0;
					vh[i, j] = 0;
					cg_mask[i, j] = false;
				}
				else
				{
					cg_mask[i, j] = true;
					b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
				}
			}
		}
		//TODO: Solve the Poisson equation to obtain vh (virtual height).
		Conjugate_Gradient(cg_mask, b, vh, lower_i - 1, upper_i + 1, lower_j - 1, upper_j + 1);
		
		//TODO: for block 2, calculate low_h.
		//TODO: then set up b and cg_mask for conjugate gradient.
		//TODO: Solve the Poisson equation to obtain vh (virtual height).
	
		//TODO: Diminish vh.
		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			if(cg_mask[i,j])
				vh[i,j]*=gamma;
		}
		
		//TODO: Update new_h by vh.
		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			if(i!=0)		new_h[i,j]+=(vh[i-1,j]-vh[i,j])*rate;
			if(i!=size-1)	new_h[i,j]+=(vh[i+1,j]-vh[i,j])*rate;
			if(j!=0)		new_h[i,j]+=(vh[i,j-1]-vh[i,j])*rate;
			if(j!=size-1)	new_h[i,j]+=(vh[i,j+1]-vh[i,j])*rate;
		}
		//Step 3
		//TODO: old_h <- h; h <- new_h;
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				old_h[i, j] = h[i, j];
				h[i, j] = new_h[i, j];
			}
		}
		//Step 4: Water->Block coupling.
		//More TODO here.
	}
	

	// Update is called once per frame
	void Update () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		Vector3[] X    = mesh.vertices;
		float[,] new_h = new float[size, size];
		float[,] h     = new float[size, size];

		//TODO: Load X.y into h.
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				h[i, j] = X[size * i + j].y; // 
			}
		}
		
		if (Input.GetKeyDown ("r")) 
		{
			//TODO: Add random water.
			float random = Random.Range(0.5f, 2.0f);
			float randomPosI = Random.Range(1f, size - 2);
			float randomPosJ = Random.Range(1f, size - 2);

			h[(int)randomPosI, (int)randomPosJ] += random;
			h[(int)randomPosI - 1, (int)randomPosJ] -= random / 4;
			h[(int)randomPosI + 1, (int)randomPosJ] -= random / 4;
			h[(int)randomPosI , (int)randomPosJ - 1] -= random / 4;
			h[(int)randomPosI, (int)randomPosJ + 1] -= random / 4;
		}
	
		for(int l=0; l<8; l++)
		{
			Shallow_Wave(old_h, h, new_h);
		}

		//TODO: Store h back into X.y and recalculate normal.
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				X[size * i + j].y = h[i, j]; // 
			}
		}

		mesh.vertices = X;
		mesh.RecalculateNormals();
	}
}
