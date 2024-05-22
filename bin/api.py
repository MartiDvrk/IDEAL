#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#   Copyright (C): MedAustron GmbH, ACMIT Gmbh and Medical University Vienna
#   This software is distributed under the terms
#   of the GNU Lesser General  Public Licence (LGPL)
#   See LICENSE for further details
# -----------------------------------------------------------------------------

# generic imports
import os
import configparser
import jwt

# ideal imports
import ideal_module as idm
import utils.condor_utils as cndr 
import utils.api_utils as ap 
import impl.dicom_functions as dcm
# api imports
from flask import Flask, request, jsonify, Response 
from flask_sqlalchemy import SQLAlchemy
from apiflask import APIFlask, HTTPTokenAuth, abort
from werkzeug.utils import secure_filename
from werkzeug.security import check_password_hash, generate_password_hash
from utils.api_schemas import SimulationRequest, Authentication, define_user_model

# Initialize sytem configuration once for all
sysconfig = idm.initialize_sysconfig(username = 'myqaion')
base_dir = sysconfig['IDEAL home']
input_dir = sysconfig["input dicom"]
log_dir = sysconfig['logging']
daemon_cfg = os.path.join(base_dir,'cfg/log_daemon.cfg')
sysconfig_path = os.path.join(base_dir,'cfg/system.cfg')
log_parser = configparser.ConfigParser()
log_parser.read(daemon_cfg)
ideal_history_cfg = log_parser['Paths']['cfg_log_file']
api_cfg = ap.get_api_cfg(log_parser['Paths']['api_cfg'])
commissioning_dir = sysconfig['commissioning']

app = APIFlask(__name__,title='IDEAL interface', version='1.0')
auth = HTTPTokenAuth(scheme='Bearer')

# api configuration
app.config['SECRET_KEY'] = os.urandom(24)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///' + api_cfg['server']['credentials db']# os.path.join(base_dir, 'database.db')
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
host_IP = api_cfg['server']['IP host']

# List of all active jobs. Members will be simulation objects
max_queue_size = 50
jobs_list = dict()
queue = ap.preload_status_overview(ideal_history_cfg,max_size=max_queue_size)

# register database 
db = SQLAlchemy(app)

User = define_user_model(db)  

@auth.verify_token
def verify_tocken(token): 
    if token is None:
        abort(401, message='Invalid token')
    try:
        data = jwt.decode(token, app.config['SECRET_KEY'], algorithms=['HS256'])
        current_user = User.query.filter_by(uid=data['public_id']).first()
    except:
        abort(401, message='Invalid token')
  
    if current_user is None:
       abort(401, message='Invalid token')

    return current_user
 

@app.route("/v1/auth", methods=['POST'])
@app.input(Authentication,location = 'headers')
def authentication(auth):
    username = auth.get('account_login')
    pwd = auth.get('account_pwd')
    user = User.query.filter_by(username=username).first()
    if not user:
        return abort(401, message='Could not verify user!', detail={'WWW-Authenticate': 'Basic-realm= "No user found!"'})  
    if not check_password_hash(user.password, pwd):
        return abort(401, message='Could not verify password!', detail= {'WWW-Authenticate': 'Basic-realm= "Wrong Password!"'})  
    token = jwt.encode({'public_id': user.uid}, app.config['SECRET_KEY'], 'HS256')

    return jsonify({'authToken': token, 'userRole':user.role, 'firstName':user.firstname, 'lastName':user.lastname}), 200 

@app.route("/v1/version")
@app.auth_required(auth)
def version():
    return idm.get_version()

@app.post("/v1/jobs")
@app.auth_required(auth)
@app.input(SimulationRequest, location='form_and_files')
def start_new_job(data):  
    # clean input directory
    ap.remove_directory_contents(input_dir)
    
    # get data from client
    rp_file = data['dicomRtPlan']
    rp_filename = secure_filename(rp_file.filename)
    rs_file = data['dicomStructureSet']
    rd_file = data['dicomRDose']
    ct_file = data['dicomCTs']

    files_dict = {}
    if ap.check_file_extension(rp_file.filename):
        files_dict["dicomRtPlan"] = 'wrong file extension'
    if ap.check_file_extension(rs_file.filename):
        files_dict["dicomStructureSet"] = 'wrong file extension'
    if ap.check_file_extension(rd_file.filename):
        files_dict["dicomRDose"] = 'wrong file extension'
    if ap.check_file_extension(ct_file.filename):
        files_dict["dicomCTs"] = 'wrong file extension'
    if files_dict:
        return jsonify(files_dict), 422
    
    arg_username = data['username']
    if not ap.check_username(sysconfig,arg_username):
        return jsonify({"username":"user not recognized"}), 400
    ref_checksum = data['configChecksum']
    arg_number_of_primaries_per_beam = data['numberOfParticles']
    arg_percent_uncertainty_goal = data['uncertainty']
    
    phantom = None 
    if 'phantom' in data:
        phantom = data['phantom']
    
    data_checksum = ap.sha1_directory_checksum(commissioning_dir,sysconfig_path)
    if data_checksum != ref_checksum:
        return jsonify({"configChecksum":"Configuration has changed from frozen original one"}), 503
    
    datadir, rp = ap.generate_input_folder(input_dir,rp_filename,arg_username)
    app.config['UPLOAD_FOLDER'] = datadir
    
    #save files in folder and  unzip dicom data
    rp_file.save(os.path.join(datadir,secure_filename(rp_file.filename)))
    rp = ap.unzip_file(datadir,rp_file.filename)
    rs_file.save(os.path.join(datadir,secure_filename(rs_file.filename)))
    rs = ap.unzip_file(datadir,rs_file.filename)
    ct_file.save(os.path.join(datadir,secure_filename(ct_file.filename)))
    cts = ap.unzip_file(datadir,ct_file.filename)
    rd_file.save(os.path.join(datadir,secure_filename(rd_file.filename)))
    rds = ap.unzip_file(datadir,rd_file.filename)
    
    # check dicom
    ok, missing_keys = dcm.verify_all_dcm_keys(datadir,rp,rs,cts,rds)
    if not ok:
        #return Response(str(missing_keys), status=422, mimetype='application/json')
        return jsonify(missing_keys), 422
    
    # create simulation object
    dicom_file = os.path.join(datadir,rp[0])
    try:
        mc_simulation = idm.ideal_simulation(arg_username,dicom_file,n_particles = arg_number_of_primaries_per_beam,
                                             uncertainty=arg_percent_uncertainty_goal, phantom = phantom)
    except Exception as e:
        abort(500, message=str(e))
        
    # # check dicom files
    # ok, missing_keys = mc_simulation.verify_dicom_input_files()
    
    # Get job UID
    jobID = mc_simulation.jobId
    
    # start simulation and append to list 
    try:
        mc_simulation.start_simulation()
    except Exception as e:
        abort(500, message=str(e))
        
    # remove oldest job if the queue has reached max size
    if len(jobs_list) >= max_queue_size:
        print('POP FIRST!')
        first_job = next(iter(jobs_list))
        jobs_list.pop(first_job)
        
    jobs_list[jobID] = mc_simulation
    
    # check stopping criteria:
    try:
        mc_simulation.start_job_control_daemon()
    except Exception as e:
        abort(500, message=str(e))
        
        
    return Response(jobID, status=201, mimetype='text/plain')
    
@app.get("/v1/jobs")
@app.auth_required(auth)
def get_queue():
    for jobId in jobs_list.keys():
        cfg_settings = jobs_list[jobId].settings
        status = ap.read_ideal_job_status(cfg_settings)
        queue[jobId] = status
    return jsonify(queue)

@app.route("/v1/jobs/<jobId>", methods=['DELETE','GET'])
@app.auth_required(auth)
def stop_job(jobId):
    if jobId not in jobs_list:
        return Response('Job does not exist', status=404, mimetype='text/plain')

    if request.method == 'DELETE':
        args = request.args
        cancellation_type = args.get('cancellationType')
        # set default to soft
        if cancellation_type is None:
            cancellation_type = 'soft'
        if cancellation_type not in ['soft', 'hard']:
            return Response('CancellationType not recognized, choose amongst: soft, hard', status=400, mimetype='text/plain')
        
        cfg_settings = jobs_list[jobId].settings
        status = ap.read_ideal_job_status(cfg_settings)
        
        if status == ap.FINISHED:
            return Response('Job already finished', status=199, mimetype='text/plain')

        if cancellation_type=='soft':
            simulation = jobs_list[jobId]
            simulation.soft_stop_simulation(simulation.cfg)
            # kill job control daemon
            try:
                daemons = cndr.get_job_daemons('job_control_daemon.py')
                cndr.kill_process(daemons[simulation.workdir])
            except:
                print('Looks like daemon is not running, not possible to kill it')
            
        if cancellation_type=='hard':
            simulation = jobs_list[jobId]
            condorId = jobs_list[jobId].condor_id
            cndr.remove_condor_job(condorId)
            # kill job control daemon
            try:
                daemons = cndr.get_job_daemons('job_control_daemon.py')
                cndr.kill_process(daemons[simulation.workdir])
            except:
                print('Looks like daemon is not running, not possible to kill it')
        
        
        
        return cancellation_type
    
    if request.method == 'GET':
        # some checks
        if jobId not in jobs_list:
            return Response('Job does not exist', status=404, mimetype='text/plain')
        if not isinstance(jobId,str):
            return Response('JobId must be a string', status=400, mimetype='text/plain')
        
        # Transfer output result upon request
        cfg_settings = jobs_list[jobId].settings
        status = ap.read_ideal_job_status(cfg_settings)
        
        if status != ap.FINISHED:
            return Response('Job not finished yet', status=409, mimetype='text/plain')
        
        outputdir = jobs_list[jobId].outputdir
        r = ap.transfer_files_to_server(outputdir,api_cfg)
        if r.status_code == 200:
            return Response('The results will be sent', status=200, mimetype='text/plain')
        else:
            return Response('Failed to transfer results', status=r.status_code, mimetype='text/plain')

@app.route("/v1/jobs/<jobId>/status", methods=['GET'])
@app.auth_required(auth)
def get_status(jobId):
    if jobId not in jobs_list:
        return Response('Job does not exist', status=404, mimetype='text/plain')
    if not isinstance(jobId,str):
        return Response('JobId must be a string', status=400, mimetype='text/plain')
    
    cfg_settings = jobs_list[jobId].settings
    status = ap.read_ideal_job_status(cfg_settings)
    return jsonify({'status': status})


if __name__ == '__main__':
    
    cert = api_cfg['server']['ssl cert']
    key = api_cfg['server']['ssl key']
    context = (cert, key)
    app.run(host=host_IP,port=5000,ssl_context=context)
    

    

# vim: set et softtabstop=4 sw=4 smartindent:
