class Robot{

	constructor(a,d){
		this.a = a;
		this.d = d;
		this.l2 = Math.sqrt(a[2] * a[2] + d[3] * d[3]);;
		this.ad2 = Math.atan2(a[2], d[3]);
	}

	static matrix(x,y,z,aDeg,bDeg,cDeg){
        let a = -aDeg * Math.PI / 180;
        let b = -bDeg * Math.PI / 180;
        let c = -cDeg * Math.PI / 180;
        let ca = Math.cos(a);
        let sa = Math.sin(a);
        let cb = Math.cos(b);
        let sb = Math.sin(b);
        let cc = Math.cos(c);
        let sc = Math.sin(c);
        let tt = [];
        tt[0] = [ca * cb, sa * cc + ca * sb * sc, sa * sc - ca * sb * cc, x];
        tt[1] = [-sa * cb, ca * cc - sa * sb * sc, ca * sc + sa * sb * cc, y];
        tt[2] = [sb, -cb * sc, cb * cc, z];
        return tt;
	}// return 3x4 matrix

	static matrix_3(aDeg,bDeg,cDeg){
		let a = -aDeg * Math.PI / 180;
        let b = -bDeg * Math.PI / 180;
        let c = -cDeg * Math.PI / 180;
        let ca = Math.cos(a);
        let sa = Math.sin(a);
        let cb = Math.cos(b);
        let sb = Math.sin(b);
        let cc = Math.cos(c);
        let sc = Math.sin(c);
        let tt = [];
        tt[0] = [ca * cb, sa * cc + ca * sb * sc, sa * sc - ca * sb * cc];
        tt[1] = [-sa * cb, ca * cc - sa * sb * sc, ca * sc + sa * sb * cc];
        tt[2] = [sb, -cb * sc, cb * cc];
        return tt;
	}// return 3x3 matrix

	static mul34_m(a,b){
		let re = [[],[],[]];
		for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 4; j++) {
                let b3j = (j === 3 ? 1 : 0);                
                re[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j] + a[i][3] * b3j;
            }
        }
        return re;
	}

	static mul34_p(a,b){
		let re = [[],[],[]];
        for (let i = 0; i < 3; i++)
            re[i] = a[i][0] * b[0] + a[i][1] * b[1] + a[i][2] * b[2] + a[i][3];
        return re;
	}

	static inverse34(m){
		let v = [[],[],[]];
		v[0][0] = -m[1][2] * m[2][1] + m[1][1] * m[2][2];
        v[0][1] = m[0][2] * m[2][1] - m[0][1] * m[2][2];
        v[0][2] = -m[0][2] * m[1][1] + m[0][1] * m[1][2];
        v[0][3] = m[0][3] * m[1][2] * m[2][1] - m[0][2] * m[1][3] * m[2][1] - m[0][3] * m[1][1] * m[2][2] + m[0][1] * m[1][3] * m[2][2] + m[0][2] * m[1][1] * m[2][3] - m[0][1] * m[1][2] * m[2][3];
        v[1][0] = m[1][2] * m[2][0] - m[1][0] * m[2][2];
        v[1][1] = -m[0][2] * m[2][0] + m[0][0] * m[2][2];
        v[1][2] = m[0][2] * m[1][0] - m[0][0] * m[1][2];
        v[1][3] = m[0][2] * m[1][3] * m[2][0] - m[0][3] * m[1][2] * m[2][0] + m[0][3] * m[1][0] * m[2][2] - m[0][0] * m[1][3] * m[2][2] - m[0][2] * m[1][0] * m[2][3] + m[0][0] * m[1][2] * m[2][3];
        v[2][0] = -m[1][1] * m[2][0] + m[1][0] * m[2][1];
        v[2][1] = m[0][1] * m[2][0] - m[0][0] * m[2][1];
        v[2][2] = -m[0][1] * m[1][0] + m[0][0] * m[1][1];
        v[2][3] = m[0][3] * m[1][1] * m[2][0] - m[0][1] * m[1][3] * m[2][0] - m[0][3] * m[1][0] * m[2][1] + m[0][0] * m[1][3] * m[2][1] + m[0][1] * m[1][0] * m[2][3] - m[0][0] * m[1][1] * m[2][3];
        return v;
	}

	static ABCby3Point(_dx, _dxy){
		var dx = LA.normalize(_dx);
        var tt = LA.mul(dx, LA.dot(_dxy, dx));
        var _dy = LA.sub(_dxy, tt);
        var dy = LA.normalize(_dy);
        var dz = LA.cross(dx, dy);
        var cacb = dx[0];
        var sacb = -dx[1];
        var sb = dx[2];
        var cbsc = -dy[2];
        var cbcc = dz[2];
        var cb = Math.sqrt(1 - sb * sb);
        var a = Math.atan2(sacb, cacb) * -180 / Math.PI;
        var b = Math.atan2(sb, cb) * -180 / Math.PI;
        var c = Math.atan2(cbsc, cbcc) * -180 / Math.PI;
        return [a, b, c];
	}

	static ABC(){
		let sb = m[2][0];
        let cb = Math.sqrt(1 - sb * sb);
        let ca = m[0][0];
        let sa = -m[1][0];
        let cc = m[2][2];
        let sc = -m[2][1];
        let a = Math.atan2(sa, ca) * -180 / Math.PI;
        let b = Math.atan2(sb, cb) * -180 / Math.PI;
        let c = Math.atan2(sc, cc) * -180 / Math.PI;
        return [a, b, c];
	}

	static flipABC_all(abc){
		return this.flipABC(abc[0], abc[1], abc[2]);
	}

	static flipABC(a,b,c){
		let na = a > 0 ? (a - 180) : (a + 180);
        let nb = b > 0 ? (180 - b) : (-180 - b);
        let nc = c > 0 ? (c - 180) : (c + 180);
        return [na, nb, nc];
	}

	inverse(v, A4_deg, base, tool){
		let Pos = Robot.matrix(v[0], v[1], v[2], v[3], v[4], v[5]);
		let goal;
		if(base==null)
			goal=Pos;
		else
			goal= Robot.mul34_m(base, Pos); // in WORLD frame BASE*POS(in BASE)= GOAL
		let T6;
		if (tool == null) {
			T6 = goal;
		} else {
			let intool = Robot.inverse34(tool);
			T6 = Robot.mul34_m(goal, intool); // T6*TOOL = GOAL= BASE*POS(in BASE)
		}
		let inreDeg = this.inverse_2(T6, A4_deg);// (T6, thetaDeg[3]);
		return inreDeg;
	}

	inverse_2(T6,A4_deg){
		let theta = [];
        let center = Robot.mul34_p(T6, [0, 0, -this.d[5]]);
        theta[0] = Math.atan2(center[1], center[0]);
        let ll = Math.sqrt(center[0] * center[0] + center[1] * center[1]);
        let p1 = [this.a[0] * center[0] / ll, this.a[0] * center[1] / ll, this.d[0]];
        let l3 = LA.dist(center, p1);
        let l1 = this.a[1];
        let beta = Math.acos((l1 * l1 + l3 * l3 - this.l2 * this.l2) / (2 * l1 * l3));
        let ttl = Math.sqrt((center[0] - p1[0]) * (center[0] - p1[0]) + (center[1] - p1[1]) * (center[1] - p1[1]));
        if (p1[0] * (center[0] - p1[0]) < 0)
            ttl = -ttl;
        let al = Math.atan2(center[2] - p1[2], ttl);
        theta[1] = beta + al;
        let gama = Math.acos((l1 * l1 + this.l2 * this.l2 - l3 * l3) / (2 * l1 * this.l2));
        theta[2] = gama - this.ad2 - 0.5 * Math.PI;
        let arr = [];
        let c = [];
        let s = [];
        for (let i = 0; i < 3; i++) {
            c[i] = Math.cos(theta[i]);
            s[i] = Math.sin(theta[i]);
        }
        arr[0] = [c[0] * (c[1] * c[2] - s[1] * s[2]), s[0], c[0] * (c[1] * s[2] + s[1] * c[2]), c[0] * (this.a[2] * (c[1] * c[2] - s[1] * s[2]) + this.a[1] * c[1]) + this.a[0] * c[0]];
        arr[1] = [s[0] * (c[1] * c[2] - s[1] * s[2]), -c[0], s[0] * (c[1] * s[2] + s[1] * c[2]), s[0] * (this.a[2] * (c[1] * c[2] - s[1] * s[2]) + this.a[1] * c[1]) + this.a[0] * s[0]];
        arr[2] = [s[1] * c[2] + c[1] * s[2], 0, s[1] * s[2] - c[1] * c[2], this.a[2] * (s[1] * c[2] + c[1] * s[2]) + this.a[1] * s[1] + this.d[0]];
        let in123 = Robot.inverse34(arr);
        let mr = Robot.mul34_m(in123, T6);
        let c5 = mr[2][2];
        if (Math.abs(c5 - 1) < 1.0E-6) {
            let A4 = -Math.PI * A4_deg / 180;
            let c4 = Math.cos(A4);
            let s4 = Math.sin(A4);
            let s6 = c4 * mr[1][0] - s4 * mr[0][0];
            let c6;
            if (Math.abs(c4) > Math.abs(s4))
                c6 = (mr[0][0] + s4 * s6) / c4;
            else
                c6 = (mr[1][0] - c4 * s6) / s4;
            theta[3] = A4;
            theta[4] = 0;
            theta[5] = Math.atan2(s6, c6);
            if (Math.abs(c6) > 1 || Math.abs(s6) > 1)
                throw new Error();
        }
        else {
            let ang = Math.atan2(mr[1][2], mr[0][2]);
            theta[3] = ang;
            theta[4] = Math.acos(c5);
            theta[5] = Math.atan2(mr[2][1], -mr[2][0]);
        }
        let inreDeg = this.rad_deg(theta);
        return inreDeg;
	}

	deg_rad(ds){
		let rd = [];
        for (let i = 0; i < 6; i++)
            rd[i] = ds[i] * Math.PI / 180;
        rd[2] -= 0.5*Math.PI;
        rd[5] += Math.PI;
        for (let i = 0; i < 6; i++)
            rd[i] = -rd[i];
        return rd;
	}

	rad_deg(ds){
		let rd = [];
        for (let i = 0; i < 6; i++)
            rd[i] = -ds[i];
        rd[2] += 0.5*Math.PI;
        rd[5] -= Math.PI;
        for (let i = 0; i < 6; i++)
            rd[i] = rd[i] * 180 / Math.PI;
        return rd;
	}

	forward_1(degs){
		let ts = this.deg_rad(degs);
        let c = [];
        let s = [];
        for (let i = 0; i < 6; i++) {
            c[i] = Math.cos(ts[i]);
            s[i] = Math.sin(ts[i]);
        }
        let m123 = [];
        m123[0] = [c[0] * (c[1] * c[2] - s[1] * s[2]), s[0], c[0] * (c[1] * s[2] + s[1] * c[2]), c[0] * (this.a[2] * (c[1] * c[2] - s[1] * s[2]) + this.a[1] * c[1]) + this.a[0] * c[0]];
        m123[1] = [s[0] * (c[1] * c[2] - s[1] * s[2]), -c[0], s[0] * (c[1] * s[2] + s[1] * c[2]), s[0] * (this.a[2] * (c[1] * c[2] - s[1] * s[2]) + this.a[1] * c[1]) + this.a[0] * s[0]];
        m123[2] = [s[1] * c[2] + c[1] * s[2], 0, s[1] * s[2] - c[1] * c[2], this.a[2] * (s[1] * c[2] + c[1] * s[2]) + this.a[1] * s[1] + this.d[0]];
        let m456 = [];
        m456[0] = [c[3] * c[4] * c[5] - s[3] * s[5], -c[3] * c[4] * s[5] - s[3] * c[5], c[3] * s[4], c[3] * s[4] * this.d[5]];
        m456[1] = [s[3] * c[4] * c[5] + c[3] * s[5], -s[3] * c[4] * s[5] + c[3] * c[5], s[3] * s[4], s[3] * s[4] * this.d[5]];
        m456[2] = [-s[4] * c[5], s[4] * s[5], c[4], c[4] * this.d[5] + this.d[3]];
        let arr = this.mul34_m(m123, m456);
        return arr;
	}

	forward_3(degs,base,tool){
		let T6 = this.forward_1(degs);
        let goal;
        if (tool == null) {
            goal = T6;
        }
        else {
            goal = Robot.mul34_m(T6, tool);
        }
        let pos;
        if (base == null) {
            pos = goal;
        }
        else {
            let inbase = Robot.inverse34(base);
            pos = Robot.mul34_m(inbase, goal);
        }
        let as = Robot.ABC(pos);
        return [pos[0][3], pos[1][3], pos[2][3], as[0], as[1], as[2]];
	}

	forwardSequence(degs){
		let ts = this.deg_rad(degs);
        let c = [];
        let s = [];
        for (let i = 0; i < 6; i++) {
            c[i] = Math.cos(ts[i]);
            s[i] = Math.sin(ts[i]);
        }
        let a0 = [];
        a0[0] = [c[0], 0, s[0], this.a[0] * c[0]];
        a0[1] = [s[0], 0, -c[0], this.a[0] * s[0]];
        a0[2] = [0, 1, 0, this.d[0]];
        let a1 = [];
        a1[0] = [c[1], -s[1], 0, this.a[1] * c[1]];
        a1[1] = [s[1], c[1], 0, this.a[1] * s[1]];
        a1[2] = [0, 0, 1, 0];
        let a2 = [];
        a2[0] = [c[2], 0, s[2], this.a[2] * c[2]];
        a2[1] = [s[2], 0, -c[2], this.a[2] * s[2]];
        a2[2] = [0, 1, 0, 0];
        let a3 = [];
        a3[0] = [c[3], 0, -s[3], 0];
        a3[1] = [s[3], 0, c[3], 0];
        a3[2] = [0, -1, 0, this.d[3]];
        let a4 = [];
        a4[0] = [c[4], 0, s[4], 0];
        a4[1] = [s[4], 0, -c[4], 0];
        a4[2] = [0, 1, 0, 0];
        let a5 = [];
        a5[0] = [c[5], -s[5], 0, 0];
        a5[1] = [s[5], c[5], 0, 0];
        a5[2] = [0, 0, 1, this.d[5]];
        let M = [];
        M[0] = a0;
        M[1] = Robot.mul34_m(M[0], a1);
        M[2] = Robot.mul34_m(M[1], a2);
        M[3] = Robot.mul34_m(M[2], a3);
        M[4] = Robot.mul34_m(M[3], a4);
        M[5] = Robot.mul34_m(M[4], a5);
        return M;
	}





}