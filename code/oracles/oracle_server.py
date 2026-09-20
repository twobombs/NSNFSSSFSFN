#/usr/bin/env python3

import socket
import re

from queue import Queue
from threading import Thread
from multiprocessing.pool import ThreadPool

from oracle_helpers import *

MAX_UNACCEPTED_CONN = 5

class OracleServerException(Exception):
    pass

class OracleServer:
    def __init__(self, logging, host, port,
                 exec_threads_nr = 1,
                 max_unaccepted_connections = MAX_UNACCEPTED_CONN):
        self.host       = host
        self.port       = port
        self.logging    = logging

        self.exec_threads_nr            = exec_threads_nr
        self.max_unaccepted_connections = MAX_UNACCEPTED_CONN

        self.queue      = Queue()
        self.buf_size   = BUF_SIZE

        self.ops = {
            OP_KEYGEN: (self.keygen, self._keygen_argparse),
            OP_SIGN: (self.sign, self._sign_argparse),
            OP_VERIFY: (self.verify, self._verify_argparse),
            OP_WRAP: (self.wrap, self._wrap_argparse)
        }

    def _exec_thread(self, id):
        log_tag = f"[_exec_thread] "
        self.logging.info(f"{log_tag}Start executor thread #{id}.")

        while True:
            self.logging.info(f"{log_tag}Waiting for next item in queue.")
            item = self.queue.get()
            self.logging.info(f"{log_tag}Retrieved following item from queue: {item}")

            if item == END_OF_QUEUE:
                self.logging.info(f"{log_tag}Queue was shut down, stopping executor.")
                if hasattr(self, "last_conn"):
                    self.last_conn.sendall(build_query(OP_CLOSE_CONN, 0, 0) + b"\n")
                self.queue.put(END_OF_QUEUE)
                break

            try:
                data, conn = item
                self.last_conn = conn
                op, ident, query = parse_query(data)
                self.logging.info(f"{log_tag}Exec thread received operation {op}, id: {ident}, query: {query}")

                if op not in self.ops:
                    self.logging.error(f"{log_tag}Received query for unsupported operation {op} (id: {ident}, query: {query}, raw data: {data})")
                    continue

                op_fn, op_argparse = self.ops[op]
                success, args = op_argparse(query)
                if success:
                    ans = op_fn(*args)
                    resp = build_query(op + 100, ident, ans)
                    conn.sendall(resp + b"\n")
                    self.logging.info(f"{log_tag}Sent response {resp}")
                else:
                    self.logging.warning(f"{log_tag}Failed to parse the query {query} for operation {op}")
            except Exception as e:
                self.logging.error(f"{log_tag}Failed to execute operation {op} with error:\n{e}")

        self.logging.info(f"{log_tag}Stopped executor thread #{id}.")

    def _receiver_thread(self, id, conn):
        log_tag = f"[_receiver_thread] "
        self.logging.info(f"{log_tag}Start receiver thread #{id}.")

        data = b""
        while True:
            try:
                while b"\n" not in data:
                    data += conn.recv(self.buf_size)

                idx = data.index(b"\n")
                curr_data = data[:idx]
                data = data[idx+1:]

                self.logging.info(f"{log_tag}Oracle server received data {curr_data}")

                if not curr_data or len(curr_data) <= 0:
                    break

                op, _, _ = parse_query(curr_data)
                if op == OP_CLOSE_CONN:
                    break

                self.queue.put((curr_data, conn))
            except Exception as e:
                self.logging.info(f"{log_tag}Caught the following exception, stopping receiver thread.\n{e}")
                break
        self.logging.info(f"{log_tag}Stopped receiver thread #{id}")

    def _run_listener(self):
        log_tag = f"[_run_listener] "

        self.logging.info(f"{log_tag}Start listener thread.")
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.bind((self.host, self.port))
        self.socket.listen(self.max_unaccepted_connections)

        self.logging.info(f"{log_tag}Oracle server listening on {self.host}:{self.port}")

        threads = []
        while True:
            try:
                self.logging.info(f"{log_tag}Waiting for connection")
                conn, addr = self.socket.accept()
                self.logging.info(f"{log_tag}Oracle server accepted connection from {addr}")

                receiver_thread = Thread(target=self._receiver_thread, args=(len(threads), conn,))
                receiver_thread.start()
                threads.append(receiver_thread)
            except Exception as e:
                self.logging.info(f"{log_tag}Caught the following exception, stopping listener thread.\n{e}")
                break
            finally:
                for thread in threads:
                    thread.join()
        self.logging.info(f"{log_tag}Stopped listener thread")

    def run(self):
        log_tag = f"[OracleServer.run] "

        self.logging.info(f"{log_tag}Oracle server starting")

        p = ThreadPool(self.exec_threads_nr)
        try:
            p.map_async(self._exec_thread, list(range(self.exec_threads_nr)))
            self._run_listener()
        except KeyboardInterrupt:
            self.logging.info(f"{log_tag}Caught user interrupt, terminating oracle...")

        finally:
            try:
                self.logging.info(f"{log_tag}Waiting for {self.queue.unfinished_tasks} operations to be finished...\
                    (Press CTRL+C again to abort immediately)")
                self.queue.put(END_OF_QUEUE)
            except KeyboardInterrupt:
                self.logging.info(f"{log_tag}Emptying query queue...")
                while not self.queue.empty():
                    self.queue.get()
                self.queue.put(END_OF_QUEUE)
            finally:
                p.close()
                p.join()
                self.close()

        self.logging.info(f"{log_tag}Oracle server stopped")

    def close(self):
        self.logging.info("[OracleServer.close] Closing socket if available.")
        if hasattr(self, "socket"):
            self.socket.close()

    def _keygen_argparse(self, query):
        m = re.match(r'(\d+):(\d+)\:(.*)', query)
        groups = m.groups()
        if len(groups) != 3:
            self.logging.error(f"[_keygen_argparse] Failed to parse keygen argument {query}, expected format 'pub_exp:modulus_bits:label' but extracted groups {groups}")
            return False, None

        pub_exp, modulus_bits, label = groups
        return True, (int(pub_exp), int(modulus_bits), label)

    def keygen(self, query):
        raise NotImplementedError("Must be implemented by child class.")

    def _sign_argparse(self, query):
        m = re.match(r'([0-9a-f]+)\:(.*)', query)
        groups = m.groups()
        if len(groups) != 2:
            self.logging.error(f"[_sign_argparse] Failed to parse sign argument {query}, expected format 'hex(data):label' but extracted groups {groups}")
            return False, None

        data_hex, label = groups
        data_bin = bytes.fromhex(data_hex)
        return True, (data_bin, label)

    def sign(self, query):
        raise NotImplementedError("Must be implemented by child class.")

    def _verify_argparse(self, query):
        m = re.match(r'([0-9a-f]+):([0-9a-f]+)\:(.*)', query)
        groups = m.groups()
        if len(groups) != 3:
            self.logging.error(f"[_verify_argparse] Failed to parse verify argument {query}, expected format 'hex(data):hex(sig):label' but extracted groups {groups}")
            return False, None

        data_hex, sig_hex, label = groups
        data_bin = bytes.fromhex(data_hex)
        sig_bin = bytes.fromhex(sig_hex)
        return True, (data_bin, sig_bin, label)

    def verify(self, query):
        raise NotImplementedError("Must be implemented by child class.")
    
    def _wrap_argparse(self, query):
        m = re.match(r'(.*)\:(.*)', query)
        groups = m.groups()
        if len(groups) != 2:
            self.logging.error(f"[_wrap_argparse] Failed to parse wrap argument {query}, expected format 'wrapping_key_label:wrapped_key_label' but extracted groups {groups}")
            return False, None

        wrapping_key_label, wrapped_key_label = groups
        return True, (wrapping_key_label, wrapped_key_label)
    
    def wrap(self, query):
        raise NotImplementedError("Must be implemented by child class.")