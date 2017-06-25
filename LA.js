class LA{
	static dist(a,b){
		let r = 0;
        for (let i = 0; i < a.length; i++)
            r += (a[i] - b[i]) * (a[i] - b[i]);
        return Math.sqrt(r);
	}

	static sub(a,b){
		let re = [];
        for (let i = 0; i < a.length; i++)
            re[i] = a[i] - b[i];
        return re;
	}

	static mag(a,b){
		let r = 0;
        for (let i = 0; i < a.length; i++)
            r += a[i] * a[i];
        r = Math.sqrt(r);
        return r;
	}

	static normalize(a){
		let r = LA.mag(a);
        let re = [];
        for (let i = 0; i < a.length; i++)
            re[i] = a[i] / r;
        return re;
	}

	static dot(a,b){
 		let re = 0;
        for (let i = 0; i < a.length; i++)
            re += a[i] * b[i];
        return re;

	}

	static mul(a,s){
		let re = [];
        for (let i = 0; i < a.length; i++)
            re[i] = a[i] * s;
        return re;
	}

	static cross(b,c){
		let re = [];
        re[0] = b[1] * c[2] - b[2] * c[1];
        re[1] = b[2] * c[0] - b[0] * c[2];
        re[2] = b[0] * c[1] - b[1] * c[0];
        return re;
	}
}