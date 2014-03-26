#include <Python.h>

double sum_array_d(int a[], int num_elements)
{
  int i;
  double sum=0;
  for (i=0; i<num_elements; i++)
    {
      sum = sum + a[i];
    }
  return(sum);
}

double _xlog2x(double d)
{
  if (d <= 0)
    return 0;
  else
    return d*log(d)/log(2.0);
}

double _choose(double n, double m)
{
  if (n<0 || m<0)
    return 0;
  else
    return exp(lgamma(n+1) - lgamma(m+1) - lgamma(n-m+1));
}

double _p(double n, double m, double p)
{
  if (n<0 || m<0 || p<0 || p>1)
    return 0;
  else
    return exp(lgamma(n+1)  - lgamma(n-m+1) - lgamma(m+1) +  m*log(p) + (n-m)*log(1.0-p));
}

double _pxlog2x(double n, double m, double p)
{
  return _p(n,m,p)*_xlog2x(m/n);
}

double _entropy(double a1[], double n1, int l1)
{
  int i;
  double n1Inv=1.0/n1;
  double h=0; 

  for (i=0;i<l1;i++){
    h-=_xlog2x(a1[i]*n1Inv);
  }

  return h;
}

double _expected_entropy(double n, double a1[], double n1, int l1)
{
  int i;
  double h=0,m=0;
  double n1inv = 1.0/n1;

  for (i=0;i<l1;i++){
    for (m=1; m<(n<a1[i] ? n : a1[i]); m++){
      h-=_pxlog2x(n,m,a1[i]*n1inv);
    }
  }

  return h;
}

double _information(double a1[], double a2[], double n1, double n2, int l)
{
  double h=_entropy(a1,n1,l);
  double havg=_expected_entropy(n1,a2,n2,l); 
  return (n1/n2)*(havg-h);
}


static PyObject* xlog2x(PyObject* self, PyObject* args)
{
  double v;
 
  if (!PyArg_ParseTuple(args, "d", &v))
    return NULL;
 
  return Py_BuildValue("d", _xlog2x(v));
}

static PyObject* choose(PyObject* self, PyObject* args)
{
  double v,w;
 
  if (!PyArg_ParseTuple(args, "dd", &v, &w))
    return NULL;
 
  return Py_BuildValue("d", _choose(v,w));
}

static PyObject* p(PyObject* self, PyObject* args)
{
  double u,v,w;
 
  if (!PyArg_ParseTuple(args, "ddd", &u, &v, &w))
    return NULL;
 
  return Py_BuildValue("d", _p(u,v,w));
}

static PyObject* pxlog2x(PyObject* self, PyObject* args)
{
  double u,v,w;
 
  if (!PyArg_ParseTuple(args, "ddd", &u, &v, &w))
    return NULL;
 
  return Py_BuildValue("d", _p(u,v,w)*_xlog2x(v/u));
}


static PyObject* entropy(PyObject* self, PyObject* args)
{
  int i,len;
  double tmp;
  PyObject* obj;
  PyObject* seq;
  PyObject* item;

  if (!PyArg_ParseTuple(args, "O", &obj))
    return NULL;

  seq = PySequence_Fast(obj, "expected a sequence");
  len = PySequence_Size(obj);

  double a1[len];
  double n1=0;

  for (i = 0; i < len; i++) {
    item = PyList_GET_ITEM(seq, i);
    if (PyFloat_Check(item)){
      tmp = PyFloat_AsDouble(item);
    }
    else if (PyInt_Check(item)){
      tmp = (double) PyInt_AsLong(item);
    }
    a1[i]=tmp;
    n1+=tmp;
  }

  Py_CLEAR(seq);
 

  return Py_BuildValue("d", _entropy(a1,n1,len)); 
}

static PyObject* expectedentropy(PyObject* self, PyObject* args)
{
  int i,len;
  double tmp, n;
  PyObject* obj;
  PyObject* seq;
  PyObject* item;

  if (!PyArg_ParseTuple(args, "Od", &obj, &n))
    return NULL;

  seq = PySequence_Fast(obj, "expected a sequence");
  len = PySequence_Size(obj);

  double a1[len];
  double n1=0;

  for (i = 0; i < len; i++) {
    item = PyList_GET_ITEM(seq, i);
    if (PyFloat_Check(item)){
      tmp = PyFloat_AsDouble(item);
    }
    else if (PyInt_Check(item)){
      tmp = (double) PyInt_AsLong(item);
    }
    a1[i]=tmp;
    n1+=tmp;
  }

  Py_CLEAR(seq);
 

  return Py_BuildValue("d",  _expected_entropy(n,a1,n1,len));
}


static PyObject* information(PyObject* self, PyObject* args)
{
  int i,len;
  double tmp=0, n1=0, n2=0;
  PyObject* obj1;
  PyObject* obj2;
  PyObject* seq1;
  PyObject* seq2;
  PyObject* item;

  if (!PyArg_ParseTuple(args, "OO", &obj1, &obj2))
    return NULL;

  seq1 = PySequence_Fast(obj1, "expected a sequence");
  seq2 = PySequence_Fast(obj2, "expected a sequence");
  len = PySequence_Size(obj1);

  double a1[len];
  double a2[len];

  for (i = 0; i < len; i++) {
    item = PyList_GET_ITEM(seq1, i);
    if (PyFloat_Check(item)){
      tmp = PyFloat_AsDouble(item);
    }
    else if (PyInt_Check(item)){
      tmp = (double) PyInt_AsLong(item);
    }
    a1[i]=tmp;
    n1+=tmp;

    item = PyList_GET_ITEM(seq2, i);
    if (PyFloat_Check(item)){
      tmp = PyFloat_AsDouble(item);
    }
    else if (PyInt_Check(item)){
      tmp = (double) PyInt_AsLong(item);
    }
    a2[i]=tmp;
    n2+=tmp;
  }

  Py_CLEAR(seq1);
  Py_CLEAR(seq2);

  return Py_BuildValue("d",  _information(a1,a2,n1,n2,len));
}

static PyMethodDef MyExtensionsMethods[] = {
  {"xlog2x", xlog2x, METH_VARARGS, "Calculate xlog2x."},
  {"choose", choose, METH_VARARGS, "Calculate choose."},
  {"p", p, METH_VARARGS, "Calculate p."},
  {"pxlog2x", pxlog2x, METH_VARARGS, "Calculate p*xlog2x."},  
  {"entropy", entropy, METH_VARARGS, "Calculate entropy of a list of (float)values."},
  {"expectedentropy", expectedentropy, METH_VARARGS, "Calculate expected entropy for n instances of a word among many."},
  {"information", information, METH_VARARGS, "Calculate information content of signal withrespect to background."},
  {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC initmyextensions(void)
{
  (void) Py_InitModule("myextensions", MyExtensionsMethods);
}
